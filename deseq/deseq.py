#!/usr/bin/env python

'''
Assumptions:
- LIMS API Token and URL are in user path
- Aggregation level cutcounts and called peak files are in LIMS

Input:

- metadata.txt file in tsv format
  agg_id, label, and condition are required.
  You can also add any additional factors you would like to include.
  as: agg_id | label | condition  | factor1 | factor2 | factor3 etc

  example:

   agg_id   label   condition   time
   5536     WT      WT          day2
   5537     DNMT3   DNMT3+      day2

    Column name for any other grouping should match experimental design:
        example:
        experimental_design = condition + time

- config.txt: Config file containing the following options.

Processing options:

- Input Master (default= None):

    Do you already have a master DHS reference, or should it be created from your input samples?
    - "none" = create from samples
    - path/to/master_dhs.bed = use this file for a reference

- Threshold (default = 0):

    Remove any rows from master DHS reference that have loci in less than X samples

- Normalization Method (default = Custom):

    DESeq2 normalization is tailored for RNAseq data and does not fit well to differential DHS loci data
    because many regions will have no DHS sites, thus making sites with DHS seem over-important.
    To deal with this issue, John Lazar has created a normalization method that uses the Geometric Mean
     or Loess smoothing methods to provide a much better fit for DHS data. This method is recommended.
    - "custom" = Use John Lazar's custom normalization method
    - "none" = Use DESeq normalization
- Geomean (default = True):

    If using custom normalization method:
    - "true" = use Geometric Mean
    - "false" = use Loess smoothing

DESeq2 options:

- Experimental design: Model used in DESeq2

    Design indicates the effects we want to measure. For example the following would examine the effect
    of condition and time on the data.
    - ex: "condition + time"
    The following would examine the effect of condition on the data, controlling for batch effect:
    - ex: "~batch + condition"
    The experimental design should only contain column names and values found in your agg_id.txt file
    - More information [here](/examples/experimental_design.R)

- Contrasts: Allow you to examine the differences in the experimental design between to groupings.

    These should be specified as a comma separated list containing (in order): "column_name,group1,group2".
    - ex: "condition,treatment,time"

-Save R objects: Save intermediate R objects for reference, and to assist in running contrasts with the
    same experiemental design more quickly. TRUE to save, FALSE to skip.

To Run:

module load bedops/2.4.20
module load bedtools/2.25.0
module load python/2.7.11
module load numpy
module load scipy/0.18.1
module load scipy/0.17.1
module load pandas
module load requests
module load atlas-lapack
module load gcc
module load R/3.2.5 # required for DEseq2
module load cairo/1.14.2 # Required for DEseq2

python deseq.py -c deseq_config.txt

'''

#################### Setup #############################

#### import packages required to load modules ####
import logging    # create log file
import os    # traverse/operations on directories and files
import subprocess as sp    # submit system calls
import time    # time/date stamp for log

### setup logging ###
'''
Setup logging here to log any module errors on loading
'''
# script log created as today_deseq2_py.log in working directory
today = time.strftime("%Y%m%d")
now = time.strftime("%Y%d%m%H%M%S")
FORMAT = "[%(name)s | Line:%(lineno)s| Function:%(funcName)s ] %(message)s"
logging.getLogger("urllib3").setLevel(logging.WARNING)
logging.getLogger("requests").setLevel(logging.WARNING)
logging.basicConfig(filename="{}_deseq2_py.log".format(today), level=logging.INFO, filemode='w',
                    format=FORMAT)
log = logging.getLogger('DESeq')

# import remaining packages
import sys
import ConfigParser as cp    # parse input config file
import argparse    #parse user arguments
import glob    # traverse directory strutctures
import pandas as pd    # parse input files
import requests    # get metadata from LIMS
import socket   # check which cluster script is being run on
import normalization.deseq_norm_v2 as norm    # John Lazar's normalization
import tests.deseq_checks as dc    # testing script for pipeline

# get path to config file
parser = argparse.ArgumentParser(description='Run DESeq pipeline for DHS expression.')
parser.add_argument('-c','--config', help='Path to deseq_config.txt',required=True)
args = parser.parse_args()

class dnaseSeq():

    # setup LIMS API URL
    try:
        if os.environ["LIMS_API_URL"]:
            api_url = os.environ["LIMS_API_URL"]
        else:
            api_url = "https://lims.altiusinstitute.org/api"
    except Exception as e:
        raise SystemExit("Make sure your LIMS API URL is setup and in your path.")

    # setup LIMS API Token
    try:
        if os.environ["LIMS_API_TOKEN"]:
            api_token = os.environ["LIMS_API_TOKEN"]
    except Exception as e:
        raise SystemExit("Make sure your LIMS API token is setup and in your path.")
    api_headers = {'Authorization': 'Token {}'.format(api_token)}

    # locate rscript path
    r_script = "{}/deseq.R".format(os.path.dirname(__file__))

    def __init__(self, agg_list, contrasts, experimental_design="condition",
                 input_master = "none", working_dir=os.getcwd(), sample_threshold=0,
                 norm_method='custom', geomean='true', input_method='peaks',
                 combine_method='FWHM', save_r_objects='true', label_col="condition", FDR=0.01):

        # set user input variables
        self.input_master = input_master # custom or pre-built master dhs list
        self.agg_list = agg_list # list of aggregation id's and metadata
        self.experimental_design = experimental_design # R experimental design parameter
        self.working_dir = working_dir # working directory for input/output
        self.sample_threshold = int(sample_threshold) # filtering will be done in R
        self.norm_method = norm_method # lazar normalization or R
        self.geomean = geomean # geomean or lowess smoothing option
        self.input_method = input_method # select between a comparison using peaks or hotpots
        self.combine_method = combine_method # count overlaps using tradition of FWHM method
        self.contrasts = contrasts # two columns to compare in samples
        self.save_r_objects = save_r_objects # boolean to save intermediate r objects
        self.label_col = label_col # name of separate column to use for labeling samples in graphs
        self.FDR = FDR # FDR cutoff for significance

        # setup output directories based on input wd
        self.input_file_dir = "{wd}/input_files/".format(wd=self.working_dir) # location of all input file symlinks
        self.count_dir = "{}/counts/".format(self.working_dir) # location of overlap counts, master dhs, normalization files
        self.temp_dir = "{}/temp".format(self.working_dir) # deleted after script is finished running
        self.results_dir = "{}/results".format(self.working_dir) # analysis results and visualizations

        # generate output file names based on wd and directory
        self.metadata_file = "{}metadata.txt".format(self.input_file_dir) # input file paths and metadata
        self.master_dhs = "{}master_dhs.bed".format(self.count_dir) # master dhs list if generated from input
        self.combined_counts = "{}/combined_counts.bed".format(self.count_dir) # combined file of overlaps between each sample and master dhs
        # self.count_order = "{count_dir}/count_order.txt".format(count_dir=self.count_dir) #
        self.norm_matrix = "{count_dir}/norm_matrix.bed".format(count_dir=self.count_dir) # binary normalization matrix between master and hotspots or peaks
        self.norm_size_factors = "{}/combined_counts.bed.size_factors".format(self.count_dir) # normalization factors generated by lazar normalization method
        self.norm_counts = "{}/combined_counts.bed.norm_counts".format(self.count_dir) # normalization counts generated by lazar normalization method

        # read in aggregate id metadata file
        self.agg_df = pd.read_csv(self.agg_list, sep="\s+", header=0) # input metadata file read in and converted to pandas df

        print("Script log created at {}/{}_deseq2_py.log \n".format(self.working_dir,today))

    @staticmethod
    def check_server():
        '''
        check if script is running on cluster or SLURM
        '''
        if socket.gethostname().startswith("sched"):
            print("Script running on SLURM cluster.")
            log.info("Script running on SLURM cluster.")
            on_slurm = "true"
        else:
            print("Script will spawn one process for each input agg_id if not run on a SLURM cluster.")
            on_slurm = "false"
        return on_slurm

    @staticmethod
    def create_outdir(output_dir):
        '''
        Check if a directory exists, and if not, create it
        :param output_dir: path to directory
        :return: directory created
        '''
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    @staticmethod
    def subprocess_cmd(command, cwd=os.getcwd(), shell=False):
        '''
        Run programs in bash via subprocess
        :param command: command string as would be run on the command line
        :param input_dir: optional directory to run command in, default cwd
        :return: runs bash command
        '''
        print ("Running \n {}".format(command))
        cmd_to_arg = command.split()
        if shell == "False":
            ps = sp.Popen(cmd_to_arg, stdout=sp.PIPE, stderr=sp.PIPE, cwd= cwd, shell=False)
        else:
            ps = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE, cwd=cwd, shell=True)
        try:
            o = ps.communicate()
            log.info(o[0])
            return o[0]
        except sp.CalledProcessError as e:
            log.error(e)
            raise SystemExit(e)

    @staticmethod
    def clean_up(dir, pattern):
        '''
        Clean up temp files and dirs
        :param dir: Directory to clean
        :param pattern: file name or wild-card style file name: ex. slurm*.out
        :return: removal of specified files in dir
        '''
        for file in glob.glob("{}/{}".format(dir, pattern)):
            os.remove(file)
            log.info("Removed temp file: {} \n".format(file))

    @staticmethod
    def check_exists(file_path):
        '''
        check that file exists and is not empty
        :param file_path:
        '''
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            print("{} exists and is not empty.".format(file_path))
        else:
            error_msg = "{} does not exist or is empty".format(file_path)
            log.error(error_msg)
            raise SystemExit(error_msg)

#################################################

    # TODO API scraping will change once LIMS endpoint is created. Please excuse the mess.
    def get_file_paths(self, agg_id):
        '''use aggregation id to get paths to cutcounts and hotspots/peaks'''
        get_counts = "{}/file/?object_id={}&object_content_type=125&purpose__slug=cutcounts-starch".format(
            self.api_url, agg_id)

        # check if input_method is using hotspots or peaks
        if self.input_method.lower() == "hotspots":
            get_input = "{}/file/?format=json&object_content_type=125&object_id={}&purpose__slug=hotspot-calls".format(
                self.api_url, agg_id)
        elif self.input_method.lower() == "peaks":
            get_input = "{}/file/?format=json&object_content_type=125&object_id={}&purpose__slug=hotspot-peaks".format(
                self.api_url, agg_id)
        else:
            log.error("Input method must be either 'hotspots' or 'peaks'.")
            raise SystemExit("Input method must be either 'hotspots' or 'peaks'.")

        try:
            v = requests.get(get_counts, headers=self.api_headers).json()
            cutcounts = v['results'][0]['path']
            p = requests.get(get_input, headers=self.api_headers).json()
            input_baseline = p['results'][0]['path']
            return cutcounts, input_baseline
        except Exception as e:
            log.error("Error found: {}. Make sure peak and cutcount files exist for AG{} and check LIMS API URL and token.".format(e, agg_id))
            raise SystemExit("Error found: {}. Make sure peak and cutcount files exist for AG{} and check LIMS API URL and token..".format(e, agg_id))

    def get_metadata(self, agg_id):
        '''
        take input agg_id's to query lims for metadata
        :param agg_id:
        :param input_file_dir:
        :return:
        '''
        try:
            get_agg = "{}/aggregation/{}".format(self.api_url, agg_id)
            r = requests.get(get_agg, headers=self.api_headers).json()
            sample_id, library = r['sample_number'], r['library_number']
            get_sample_info = "{}/sample/?number={}".format(self.api_url, r['sample_number'])
            t = requests.get(get_sample_info, headers=self.api_headers).json()
            taxonomy, assay, url = t['results'][0]['sample_taxonomy_name'], t['results'][0]['assay_name'], \
                               t['results'][0]['url']
            return sample_id, library, taxonomy, assay,url
        except Exception as e:
            log.error("Error found: {}. Make sure sample id, library number, \
                             taxonomy, assay and view url exist.".format(e))
            raise SystemExit("Error found: {}. Make sure sample id, library number, \
                             taxonomy, assay and view url exist.".format(e))

    # read in aggregation id's
    def get_api_metadata(self, agg_df):
        '''
        Use aggregation id's to query LIMS and create metadata file
        :param agg_file: text file of aggregation id's, one per line, no header
        :return: dataframe containing agg_id, sample_id, library_number, taxonomy, assay, cutcounts_path
        file metadata.txt
        '''
        meta_list = []
        # query LIMS for each aggregation id
        for agg in agg_df['agg_id']:
            if agg.lower().startswith("ag"):
                agg_id = int("{}".format(str(agg.lower()).replace('ag','')))
            cutcounts, input_baseline = self.get_file_paths(agg_id)
            meta_row = [cutcounts, input_baseline]
            meta_list.append(meta_row)
        # create dataframe of metadata
        if self.input_method.lower() == "hotspots":
            meta_cols = pd.DataFrame(meta_list, columns=['cutcounts', 'hotspots'])
            meta_df = pd.concat([agg_df, meta_cols], axis=1, join='inner')
        elif self.input_method.lower() == "peaks":
            meta_cols = pd.DataFrame(meta_list, columns=['cutcounts', 'peaks'])
            meta_df = pd.concat([agg_df, meta_cols], axis=1, join='inner')
        else:
            log.error("Input method must be either 'hotspots' or 'peaks'.")
            raise SystemExit("Input method must be either 'hotspots' or 'peaks'.")

        return meta_df

    def get_manual_input(self, agg_df):
        '''
        Check for non-LIMS file paths and API calls and gather metadata
        :param agg_df: self.agg_df pandas dataframe of input metadata text with agg id and/or filepaths
        :return: input_files/metadata.txt and dataframe to be fed to function for creating symlinks to input files
        '''

        if self.input_method.lower() == "hotspots":
            # if the metadata contains a peak or cutcounts file
            if agg_df['cutcounts'].notnull().values.any() or agg_df['hotspots'].notnull().values.any():
                # subset rows containing peak and cutcount file paths
                manual_df = agg_df[agg_df.cutcounts.notnull()]
                # check if there are any api calls
                api_df = agg_df[agg_df.cutcounts.isnull()]
                api_df = api_df.drop(['cutcounts', 'hotspots'], axis=1)

                # check that both peak and cutcount paths are provided
                if manual_df['cutcounts'].isnull().values.any() or manual_df['hotspots'].isnull().values.any():
                    raise SystemExit("Both hotspot and cutcount file paths must be provided for samples with no agg id.")
                # if both peak and cutcount files provided, check for api calls and return metadata
                else:
                    # if there are also api calls, get metadata and combine with manual metadata
                    if not api_df.empty:
                        # get metadata from LIMS for rows with agg id's
                        api_metadata = self.get_api_metadata(api_df)
                        # merge the two data frames to return all metadata
                        meta_df =  api_metadata.append(manual_df, ignore_index=True)
                        meta_df.sort_values('condition', inplace=True, ascending=True)

                        # check if input_files dir exists, and create if not
                        self.create_outdir(self.input_file_dir)
                        # write metadata to input_files/metadata.txt
                        meta_df.to_csv(self.metadata_file, sep='\t', index=False)
                        return meta_df
                    # if there are no api calls, just return manual metadata
                    else:
                        meta_df = manual_df.sort_values('condition', ignore_index=True, ascending=True)
                        # check if input_files dir exists, and create if not
                        self.create_outdir(self.input_file_dir)
                        # write metadata to input_files/metadata.txt
                        meta_df.to_csv(self.metadata_file, sep='\t', index=False)
                        return meta_df
            elif self.input_method.lower() == "peaks":
                # if the metadata contains a peak or cutcounts file
                if agg_df['cutcounts'].notnull().values.any() or agg_df['peaks'].notnull().values.any():
                    # subset rows containing peak and cutcount file paths
                    manual_df = agg_df[agg_df.cutcounts.notnull()]
                    # check if there are any api calls
                    api_df = agg_df[agg_df.cutcounts.isnull()]
                    api_df = api_df.drop(['cutcounts', 'peaks'], axis=1)

                    # check that both peak and cutcount paths are provided
                    if manual_df['cutcounts'].isnull().values.any() or manual_df['peaks'].isnull().values.any():
                        raise SystemExit(
                            "Both peak and cutcount file paths must be provided for samples with no agg id.")
                    # if both peak and cutcount files provided, check for api calls and return metadata
                    else:
                        # if there are also api calls, get metadata and combine with manual metadata
                        if not api_df.empty:
                            # get metadata from LIMS for rows with agg id's
                            api_metadata = self.get_api_metadata(api_df)
                            # merge the two data frames to return all metadata
                            meta_df = api_metadata.append(manual_df, ignore_index=True)
                            meta_df.sort_values('condition', inplace=True, ascending=True)

                            # check if input_files dir exists, and create if not
                            self.create_outdir(self.input_file_dir)
                            # write metadata to input_files/metadata.txt
                            meta_df.to_csv(self.metadata_file, sep='\t', index=False)
                            return meta_df
                        # if there are no api calls, just return manual metadata
                        else:
                            meta_df = manual_df.sort_values('condition', ignore_index=True, ascending=True)
                            # check if input_files dir exists, and create if not
                            self.create_outdir(self.input_file_dir)
                            # write metadata to input_files/metadata.txt
                            meta_df.to_csv(self.metadata_file, sep='\t', index=False)
                            return meta_df
            else:
                log.error("Input method must be either 'hotspots' or 'peaks'.")
                raise SystemExit("Input method must be either 'hotspots' or 'peaks'.")


    def collect_metadata(self, agg_df):
        # if manual file paths found
        if 'cutcounts' in agg_df.columns or 'hotspots' in agg_df.columns or 'peaks' in agg_df.columns:
            meta_df = self.get_manual_input(agg_df)

            # check that metadata file was created
            self.check_exists(self.metadata_file)

            return meta_df

        # otherwise just do api call
        else:
            meta_df= self.get_api_metadata(self.agg_df)
            print(meta_df)
            # check if input_files dir exists, and create if not
            self.create_outdir(self.input_file_dir)
            # write metadata to input_files/metadata.txt
            meta_df.to_csv(self.metadata_file, sep='\t', index=False)

            # check that metadata file was created
            self.check_exists(self.metadata_file)

            return meta_df


#################################################################################################

###########################
### Setup directories  ####
###########################

    def create_dirs(self):
        '''
        Create output directories if they do not exist
        :return: output directories created from input wd user arg
        '''
        self.create_outdir(self.temp_dir)
        self.create_outdir(self.count_dir)
        self.create_outdir(self.input_file_dir)
        self.create_outdir(self.results_dir)
        log.info("Input File Metadata: {} \n \
                    Temp Files: {} \n \
                    Counts: {} \n \
                    Results: {}". \
                 format(self.input_file_dir, self.temp_dir, self.count_dir, self.results_dir))

    def get_agg_files(self, input_df, input_col, input_file_dir):
        '''
        read in metadata file and create symlinks to input files
        symlinks generated in wd/input_files
        :param input_df:
        :param input_col:
        :param input_file_dir:
        :return: wd/input_files/metadata.txt with agg_id, library_number, file paths, assay
        '''
        log.info("Setting up output directories in {}".format(self.working_dir))
        for index, row in input_df.iterrows():
            src_file = row["{}".format(input_col)]
            dest_file = "{starch_dir}/{agg_id}_{col_name}.bed.starch"\
                           .format(starch_dir=input_file_dir, agg_id=row['agg_id'],
                                   col_name=input_col)
            if os.path.isfile(dest_file):
                log.info("Symlink already exists for: {}".format(dest_file))
                print("Symlink already exists for: {}".format(dest_file))
                pass
            else:
                try:
                    if os.path.isfile(src_file):
                        os.symlink(src_file, dest_file)
                        log.info("Symlink to {} created at {}".format(src_file, dest_file))
                    else:
                        log.error("{} does not exist".format(src_file))
                        raise SystemExit("{} does not exist".format(src_file))
                except OSError, e:
                    if e.errno != 17:
                        log.error("Error occured: {}".format(e))
                        raise SystemExit(e)
                    else:
                        time.sleep(5)
                        pass

##################################
#### Create master DHS List  #####
##################################

    def bedops_union_input(self, temp_dir, file_list):
        '''
        Union all input peak files into one file using bedops -u
        removing any mitochondrial reads
        :param temp_dir: temporary directory to create files in
        :param file_list: list of hotspot or peak files to merge, space separated
        :return: wd/temp_dir/tmp.bed file of unioned hotspot or peaks for all input aggregate peak files
        '''
        union_cmd = r"""bedops -u {file_list} | awk '{{if ($1 != "chrM") print $0;}}' > {temp_dir}/tmp.bed"""\
                .format(file_list=file_list, temp_dir=temp_dir)
        self.subprocess_cmd(union_cmd, shell=True)
        log.info("Bedops -u performed on {}. Temp files output to {}".format(file_list, temp_dir))

    def merge_union(self, temp_dir):
        '''
        Condense the unioned hotspots or peaks into merged intervals
        to avoid merging regions that are simply adjacent,
        but not overlapping
        :param temp_dir: path to temp_dir where tmp.bed from bedops_union_input()
        :return: tmpm.bed temp file with merged unions condensed and checked for overlaps vs adjacency
        '''
        merge_union_cmd = "bedops -m --range 0:-1 {temp_dir}/tmp.bed |\
                           bedops -u --range 0:1 - > {temp_dir}/tmpm.bed".format(temp_dir=temp_dir)
        self.subprocess_cmd(merge_union_cmd, shell=True)
        log.info("Bedops union files merged for overlapping regions.")

    def get_max(self, temp_dir, count_dir, master_dhs_out):
        '''
        For any ties, get the region with the highest score
        Score most likely represents read depth, but cannot find official documentation
        :param temp_dir: path to wd/temp
        :param count_dir: path to wd/counts
        :return: master_dhs.bed file of unioned hotspots or peaks with max score in overlap regions
        '''
        get_max_cmd = """
        bedmap --echo --max {temp_dir}/tmpm.bed {temp_dir}/tmp.bed \
        | bedmap --echo --echo-map {temp_dir}/tmp.bed - \
        | awk -F'|' '{{print $3"|"$1}}' \
        | awk -F'[|\\t]' '$1==$6' \
        | awk -F'|' '{{print $2}}' \
        | bedmap --echo-map --prec 0 {temp_dir}/tmpm.bed - \
        | awk -F';' '{{print $1}}' | cut -f1,2,3
        """.format(temp_dir=temp_dir, count_dir=count_dir,
                   master_dhs_out=master_dhs_out)
        with open(master_dhs_out, "w") as outfile:
            sp.call(get_max_cmd, stdout=outfile, shell=True)
            log.info("Master DHS file created at {}".format(master_dhs_out))

    def create_masterDHS(self):
        '''
        Function to run creation of master DHS using bedops_union(),
        merge_union() and get_max() functions above
        :return: master_dhs.bed
        '''
        print("Creating master dhs list.")
        if self.input_method.lower() == "hotspots":
            files = glob.glob("{}/*_hotspots.bed.starch".format(self.input_file_dir))
        elif self.input_method.lower() == "peaks":
            files = glob.glob("{}/*_peaks.bed.starch".format(self.input_file_dir))
        file_list = ' '.join(files)
        self.bedops_union_input(temp_dir=self.temp_dir, file_list=file_list)
        self.merge_union(temp_dir=self.temp_dir)
        self.get_max(temp_dir=self.temp_dir,
                     count_dir=self.count_dir, master_dhs_out=self.master_dhs)
        master_length = int(dc.check_output_length(self.master_dhs).split(' ')[0])
        if master_length > 0:
            log.info("{} lines in master_dhs.bed".format(master_length))
        else:
            log.error("No output found in master_dhs.bed")
            raise SystemExit("No output found in master_dhs.bed")

    def symlink_master(self, input_master, master_dhs):
        '''
        Create symlink to master DHS file in input_files dir
        :param master_dhs_file: user input arg for pre-existing master_dhs file
        :param master_dhs: path to wd/counts/master_dhs.bed
        :return:
        '''
        print("Creating symlink to master dhs list.")
        symlink_master_cmd = "ln -s {input_master} {master_dhs}".format(\
            input_master=input_master, master_dhs=master_dhs)
        self.subprocess_cmd(symlink_master_cmd, shell=True)

    def check_master(self):
        '''
        Check if master DHS file input by user, create if not
        :return: checks if user input master_dhs path, creates if not
        '''
        if (self.input_master.lower() == 'none') and (self.input_master is not None):
            self.create_masterDHS()
            self.check_exists(self.master_dhs)
        else:
            print ("Creating symlink to {}".format(self.master_dhs))
            self.symlink_master(self.input_master, self.master_dhs)
            self.check_exists(self.master_dhs)

###################################
####### Count Overlaps ############
###################################

    def run_slurm_jobs(self):
        '''
        Create command to count overlaps between master_dhs.bed
        and cutcounts.bed.starch on SLURM cluster
        using bedmap --count --sum
       :return: submits one job for each input file, returns list of job id's
        '''
        # check for master_dhs file
        if self.master_dhs:

            job_list = []

            # for each input agg_id in self.agg_df
            for index, row in self.agg_df.iterrows():

                agg_id = row['agg_id']
                job_name = "{agg_id}ct".format(agg_id=agg_id)
                file = "{input_dir}/{agg_id}_cutcounts.bed.starch".format(input_dir=self.input_file_dir,
                                                                          agg_id=agg_id)

                # create count overlaps command to be run on SLRUM
                slurm_cmd = r"""sbatch --job-name={job_name} --partition=queue1 \
                        --wrap="bedmap  --mem-per-cpu=10GB --prec 0 --delim '\t' \
                        --count --sum {master} {cut_file} |
                        awk '{{  if(0==\$1){{print "0"}} else{{print \$2}}  }}' > {count_dir}/{agg_id}_overlaps.bed" """\
                        .format(job_name=job_name, input_dir=self.input_file_dir,
                                master=self.master_dhs, cut_file=file,
                                count_dir=self.count_dir, agg_id=agg_id)
                p = self.subprocess_cmd(slurm_cmd, shell=True)
                job_id = p.split(" ")[3]
                job_list.append(job_id)
            return job_list

        else:
            log.info("No master dhs list found. Make sure to run check_master() function.")
            raise SystemExit("No master dhs list found. Make sure to run check_master()")

    def check_job_status(self, job_list):

        # check that job completed successfully
        complete_cmd = "sacct -P --delimiter='_' -n --format JobID,State -j {} | grep -v .bat ".format(job_list)

        # test that jobs completed successfully
        t = self.subprocess_cmd(complete_cmd, shell=True).strip().split('\n')
        for x in t:
            if x.split("_")[1] != "COMPLETED":
                raise SystemExit("Job {} had status {}.".format(x.split("_")[0], x.split("_")[1]))

    def run_jobs(self):
        '''
        Checks list of input job id's and waits until all jobs have left squeue
        :return: Pauses script execution until all running jobs have left squeue
        '''
        # create string list of jobs with no spaces for squeue command
        jobs = ','.join(self.run_slurm_jobs()).replace("\n", "")

        try:
            job_check_cmd = "squeue -h -j {}".format(jobs)
            job_check = self.subprocess_cmd(job_check_cmd, shell=True)
            while len(job_check) > 0:
                job_check = self.subprocess_cmd(job_check_cmd, shell=True)
                print(job_check)
                time.sleep(5)
            self.check_job_status(jobs)
        except Exception as e:
            raise SystemExit(e)

    def cluster_count_overlaps(self, input_dir, master_dhs_file, count_dir):

        # create list of count overlap commands
        command_list = []

        # check for master_dhs file
        if master_dhs_file:

            # for each input agg_id in self.agg_df
            for index, row in self.agg_df.iterrows():
                # build input cutcounts file names
                file = "{input_dir}/{agg_id}_cutcounts.bed.starch".format(input_dir=input_dir,
                                                                          agg_id=row['agg_id'])

                # build count overlaps command
                server_cmd = r""" nice --adjustment=5 \
                            bedmap --prec 0 --delim '\t' \
                            --count --sum {master} {cut_file} |
                            awk '{{  if(0==$1){{print "0"}} else{{print $2}}  }}' > \
                            {count_dir}/{agg_id}_overlaps.bed
                            """.format(master=self.master_dhs, cut_file=file,
                                       count_dir=count_dir, agg_id=row['agg_id'])

                # add command to command_list
                command_list.append(server_cmd)

            return command_list

        else:
            log.info("No master dhs list found. Make sure to run check_master() function.")
            raise SystemExit("No master dhs list found. Make sure to run check_master()")

    def run_count_overlaps(self, on_slurm):
        # if script is being run on slurm use sbatch wrap
        if on_slurm == 'true':

            # for each input agg_id in self.agg_df
            self.run_jobs()

        else:
            try:
                overlap_cmds = self.cluster_count_overlaps(input_dir=self.input_file_dir,
                                                           master_dhs_file=self.master_dhs,
                                                           count_dir=self.count_dir)
                # run in parallel
                processes = [sp.Popen(cmd, shell=True,
                                      stdout=sp.PIPE,
                                      stderr=sp.STDOUT, ) for cmd in overlap_cmds]
                # wait for completion
                for p in processes: p.wait()
                process_output, _ = p.communicate()
                log.info(process_output)
            except Exception as e:
                log.info('Exception occured: {}'.format(e))
                raise SystemExit('Exception occured: {}'.format(e))

#########################################
####### Join Overlap Counts ############
#########################################

    def create_merger_list(self, count_dir):
        '''
        Create list of overlap count files in same order as agg_id.txt
        :param count_dir: path to overlaps.bed files
        :return: in-memory list of overlaps.bed files to merge
        '''
        merge_list = []
        for index, row in self.agg_df.iterrows():
            # build cutcounts file name
            file = "{count_dir}/{agg_id}_overlaps.bed".format(count_dir=count_dir,
                                                                      agg_id=row['agg_id'])
            merge_list.append(file)
        return merge_list


    def merge_beds(self, count_dir, master_dhs):
        '''
        Create merged file containing one column per input file,
        with annotations from master_dhs list
        :param bed_dir: path to overlaps.bed files
        :return: annotated_counts.bed
        '''
        merger_list = self.create_merger_list(count_dir)
        merged_file_list = ' '.join(merger_list)
        print(merger_list)
        merge_cmd = r"""paste -d'\t' {master} {files} > {count_dir}/annotated_counts.bed""".format(
            master=master_dhs, files=merged_file_list, count_dir=count_dir)
        self.subprocess_cmd(merge_cmd, shell=True)
        remove_annotations_cmd = "cat {}/annotated_counts.bed | cut -f4- > {}".format(self.count_dir, self.combined_counts)
        self.subprocess_cmd(remove_annotations_cmd, shell=True)
        log.info("Bed files merged: \n {}".format(merge_cmd))
        annotated_length = int(dc.check_output_length("{}annotated_counts.bed".format(count_dir)).split(' ')[0])
        if annotated_length > 0:
            log.info("{} lines in annotated_counts.bed".format(annotated_length))
        else:
            log.error("No data found in annotated_counts.bed")
            raise SystemExit("No data found in annotated_counts.bed")

##################################
####### Filter Counts ############
##################################

    # def filter_merged_counts(self, count_dir):
    #     '''
    #     Filter out rows with dhs not found in at least 1 sample
    #     :param input_master: path to master_dhs.bed
    #     :param combined_counts_file: path to annotated_counts.bed
    #     :param sample_threshold: minimum number of samples dhs must be found in
    #     :return: filtered and unfiltered master_dhs and annotated counts
    #     '''
    #
    #     # read in annotated combined counts file
    #     annotation_file = "{}annotated_counts.bed".format(count_dir)
    #     annotated_counts = pd.read_csv(annotation_file, sep="\t", header=None)
    #
    #     log.info("/n ############## Filtering Counts ################## /n")
    #
    #     # make directory for unfiltered output
    #     self.create_outdir("{}/unfiltered_counts/".format(count_dir))
    #     # copy unfiltered files to unfiltered_counts as backup
    #     copy_unfiltered_cmd = "cp {out_dir}/annotated_counts.bed \
    #                 {out_dir}/unfiltered_counts/unfiltered_annotated_counts.bed; \
    #                 cp {out_dir}/master_dhs.bed {out_dir}/unfiltered_counts/unfiltered_masterdhs.bed".format(
    #                 out_dir=count_dir)
    #     self.subprocess_cmd(copy_unfiltered_cmd, shell=True)
    #
    #     # filter the annotated data frame according to input sample_threshold
    #     filtered_df = annotated_counts[(annotated_counts.ix[:,3:] > 0).sum(axis=1) >= self.sample_threshold]
    #
    #
    #     # save filtered counts as filtered_annotated_counts.bed
    #     filtered_df.to_csv("{}/annotated_counts.bed".format(count_dir),
    #                        sep="\t", header=None, index=None)
    #
    #     # save the filtered counts without annotation for normalization
    #     filtered_counts = filtered_df.ix[:,3:]
    #     filtered_counts.to_csv("{}/combined_counts.bed".format(count_dir),
    #                        sep="\t", header=None, index=None)
    #
    #     # subset the filtered dataframe to create a filtered master_dhs_list
    #     filtered_master = filtered_df.ix[:,0:2]
    #     filtered_master.to_csv("{}/master_dhs.bed".format(count_dir), sep="\t", header=None, index=None)
    #     log.info("Filtered counts saved to: {out_dir} \n Original files saved to: \
    #              {out_dir}/unfiltered_counts".format(out_dir=count_dir))
    #
    #     filtered_length = int(dc.check_output_length("{}combined_counts.bed".format(count_dir)).split(' ')[0])
    #     if filtered_length > 0:
    #         log.info("{} lines in combined_counts.bed".format(filtered_length))
    #     else:
    #         log.error("No data found in combined_counts.bed")
    #         raise SystemExit("No data found in combined_counts.bed")

#####################################
####### Normalize Counts ############
#####################################

    def create_norm_list(self):
        '''
        Create list of overlap count files in same order as agg_id.txt
        :param count_dir: path to overlaps.bed files
        :return: in-memory list of overlaps.bed files to merge
        '''
        merge_list = []
        for index, row in self.agg_df.iterrows():
            # build cutcounts file name
            if self.input_method.lower() == "hotspots":
                file = "{input_dir}/{agg_id}_hotspots.bed.starch".format(input_dir=self.input_file_dir,
                                                                      agg_id=row['agg_id'])
            elif self.input_method.lower() == "peaks":
                file = "{input_dir}/{agg_id}_peaks.bed.starch".format(input_dir=self.input_file_dir,
                                                                      agg_id=row['agg_id'])
            merge_list.append(file)
        return merge_list

    def create_normalization_matrix(self, out_dir, master_dhs_list):
        '''
        create binary matrix used for John Lazar's normalizaiton script
        by reading in the annotated counts file and converting to binary
        :param input_dir: path to input peak files
        :param out_dir: path to wd/counts
        :param master_dhs_list: path to master_dhs.bed
        :return: norm_matrix.bed
        '''

        merger_list = self.create_norm_list()
        merged_file_list = ' '.join(merger_list)
        matrix_cmd = "i=0; for file in {merged_list}; do ((i++)) ; if [ \"$i\" == \"1\" ]; \
                then bedmap --indicator {master} $file > {count_dir}norm_matrix.bed ; \
                else bedmap --indicator {master} $file \
                | paste {count_dir}norm_matrix.bed - > temp.txt\
                ; mv temp.txt {count_dir}norm_matrix.bed ; fi ; done".format(
            merged_list=merged_file_list, master=master_dhs_list, count_dir=out_dir)
        self.subprocess_cmd(matrix_cmd, shell=True)

    def run_normalization(self):
        '''
        If custom normalization method chosen,
        run John Lazar's normalization script
        :return: wd/counts/combined_counts.bed.norm_counts,
            wd/counts/combined_counts.bed.size_factors
        '''
        if self.norm_method.lower() == 'custom':
            log.info("\n ############## Normalization   ################## \n")
            self.create_normalization_matrix(out_dir=self.count_dir,
                                             master_dhs_list=self.master_dhs)
            norm.run_normalization_file(self.combined_counts, self.norm_matrix,
                                        self.geomean, self.combined_counts)
            log.info("Normalization method: John Lazar's normalization script, Geomean={}".format(self.geomean))
            dc.check_norm_matrix(self.norm_matrix)
            normalized_length = int(dc.check_output_length("{}combined_counts.bed.size_factors".format(self.count_dir)).split(' ')[0])
            if normalized_length > 0:
                log.info("{} lines in combined_counts.bed.size_factors".format(normalized_length))
            else:
                log.error("No data found in combined_counts.bed.size_factors")
                raise SystemExit("No data found in combined_counts.bed.size_factors")

        else:
            log.info("Normalization not performed.")
            pass

########################################
####### Clean up temp files ############
########################################

    def deseq_clean(self):
        print("Removing temp files and /temp directory")
        self.clean_up("{}".format(self.working_dir), "slurm*.out")
        rm_temp_dir = "rm -rf {}".format(self.temp_dir)
        self.subprocess_cmd(rm_temp_dir, shell=True)
        log.info("Removed temp files in {}".format(self.temp_dir))

###############################
####### Run DESeq2 ############
###############################

    def run_deseq2(self):
        results_dir = "{}/results/".format(working_dir)
        self.create_outdir(results_dir)
        r_cmd = "Rscript {rscript} -d \"{design}\" -m {master_dhs_list} -t {threshold}\
        -c {combined_counts} -i {metadata_file} -f {size_factors} -o {results_dir} -n {input_method}\
        -a {contrasts} -s {save_rds} -l {label_col} -p {fdr}".format(rscript=self.r_script,
                                    design=self.experimental_design, master_dhs_list=self.master_dhs,
                                    threshold=self.sample_threshold, combined_counts=self.combined_counts,
                                    metadata_file=self.agg_list, size_factors=self.norm_size_factors,
                                    results_dir=results_dir, input_method= self.input_method, contrasts=self.contrasts,
                                    save_rds=self.save_r_objects, label_col=self.label_col, fdr=self.FDR)

        print("{}".format(r_cmd))

        p= sp.Popen(r_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
        log.info("Running DESeq2")
        stdout, stderr = p.communicate()

        if stdout:
            log.info("Running DESeq2 with parameters: \n{}".format(r_cmd))

        if stderr:
            log.error(stderr)
            if p.returncode != 0:
                raise SystemExit("System exited with return code: {}. See log for details.".format(p.returncode))
        else:
            print("DHS analysis successfully completed.")

###############################
####### Run script ############
###############################

    def main(self):

        on_slurm = self.check_server()

        log.info("\n ############## Verifying input user args   ################## \n")
        print("Verifying input user args... \n")
        # Check that input args are correct using deseq_checks.py
        dc.main(master_dhs_list=self.input_master, experimental_design=self.experimental_design,
                working_dir= self.working_dir, sample_threshold=self.sample_threshold,
                norm_method=self.norm_method, geomean=self.geomean, contrasts=self.contrasts,
                label_col=self.label_col, agg_df=self.agg_df)

        # create output directories if they do not exist
        log.info("\n ############## Creating output directories   ################## \n")
        print("Creating output directories... \n")
        self.create_dirs()

        # collect metadata from agg_id file and symlink
        log.info("\n ############## Gathering metadata   ################## \n")
        print("Collecting metadata from input aggregate id's... \n")
        meta_df = self.collect_metadata(self.agg_df)
        if self.input_method.lower() == "hotspots":
            self.get_agg_files(meta_df, 'cutcounts', self.input_file_dir)
            self.get_agg_files(meta_df, 'hotspots', self.input_file_dir)
        elif self.input_method.lower() == "peaks":
            self.get_agg_files(meta_df, 'cutcounts', self.input_file_dir)
            self.get_agg_files(meta_df, 'peaks', self.input_file_dir)
        else:
            log.error("input_method must be either peaks or hotspots. \n")
            raise SystemExit("input_method must be either peaks or hotspots.")
        dc.check_symlinks(agg_df=self.agg_df, wd=self.working_dir)

        # check for input master dhs list, and create from input samples if not exists
        log.info("\n ############## Master DHS List   ################## \n")
        print("Checking master dhs list... \n")
        self.check_master()

        # count overlaps between each input aggregation and master dhs list
        log.info("\n ############## Count Overlaps   ################## \n")
        print("Counting overlaps against master dhs list... \n")
        self.run_count_overlaps(on_slurm=on_slurm)

        # merge overlap counts files
        log.info("\n ############## Merge Overlaps   ################## \n")
        print("Merging overlap counts files...\n")
        self.merge_beds(self.count_dir, self.master_dhs)

        # if normalization = custom, create binary matrix and normalize with normalization_script.py
        print("Checking for input normalization options... \n")
        self.run_normalization()

        # clean up temp files
        log.info("############## Cleanup   ##################")
        print("Cleaning up temp files... \n")
        self.deseq_clean()

        # run R script run_dnaseq.R with input params
        log.info("\n ############## Running DESeq  ################## \n")
        print("Running DESeq2... \n")
        self.run_deseq2()

if __name__ == "__main__":

    # parse config file with default options
    config = cp.RawConfigParser({"working_dir":os.getcwd(),
                                 "input_master":"none",
                                 "sample_threshold":0,
                                 "norm_method":"custom",
                                 "geomean":"true",
                                 "input_method":"peaks",
                                 "save_r_objects":"true",
                                 "label_col":"condition",
                                 "FDR":0.01
                                 })
    configFilePath = args.config
    config.read(configFilePath)

    # set input args
    working_dir = config.get("paths", "working_dir")
    input_master_dhs = config.get("paths", "input_master")
    aggregations_list = config.get("paths", "agg_list")
    threshold = config.get("transforms", "sample_threshold")
    norm_method = config.get("transforms", "norm_method")
    geomean = config.get("transforms", "geomean")
    input_method = config.get("transforms", "input_method")
    FDR = config.get("transforms", "fdr")
    exp_design = config.get("deseq_options", "experimental_design")
    contrasts = config.get("deseq_options", "contrasts")
    save_r_objects = config.get("deseq_options", "save_r_objects")
    label_col = config.get("deseq_options", "label_col")

    # instantiate class with input args
    ds = dnaseSeq(input_master=input_master_dhs, agg_list=aggregations_list,
                  experimental_design=exp_design, working_dir=working_dir, sample_threshold=threshold,
                  norm_method=norm_method, geomean=geomean, input_method=input_method, contrasts=contrasts,
                  save_r_objects=save_r_objects, label_col=label_col, FDR=FDR)

    # run main method
    ds.main()
