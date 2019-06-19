#!/usr/bin/env python

import logging
import os
import pandas as pd
import subprocess as sp
import time

# script log created as today_deseq2_py.log in working directory
today = time.strftime("%Y%d%m")
FORMAT = "[%(name)s | Line:%(lineno)s| Function:%(funcName)s ] %(message)s"
logging.getLogger("urllib3").setLevel(logging.WARNING)
logging.getLogger("requests").setLevel(logging.WARNING)
logging.basicConfig(filename="{}_deseq2_py.log".format(today), level=logging.INFO, filemode='w',
                    format=FORMAT)
log = logging.getLogger('DESeq_checks')

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
        ps = sp.Popen(cmd_to_arg, stdout=sp.PIPE, stderr=sp.PIPE, cwd=cwd, shell=shell)
    else:
        ps = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE, cwd=cwd, shell=shell)
    try:
        o = ps.communicate()
        log.info(o)
        return (o[0])
    except sp.CalledProcessError as e:
        log.info(e)
        raise SystemExit(e)

def compare_file_length_dir(input_dir, input_extension):
    compare_cmd = "for file in {}*{}; do wc -l $file| cut -f1 -d' '; done".format(input_dir, input_extension)
    results = subprocess_cmd(compare_cmd, shell=True)
    return results

def check_output_length(input_file_path):
    length_cmd = "wc -l {}".format(input_file_path)
    results = subprocess_cmd(length_cmd, shell=True)
    return results

def check_norm_matrix(binary_matrix):
    non_zero_count = subprocess_cmd("cat {} | grep -v 0| wc -l".format(binary_matrix), shell=True)
    if non_zero_count[0] == 0:
        log.error("Normalization matrix contains only 0's.")
        raise SystemExit("Normalization matrix containts only 0's.")
    else:
        log.info("Normalization matrix: passed.")

def check_contrasts(input_contrast_string, agg_df):

    # format input contrasts
    contrasts = [x.replace('"', '') for x in input_contrast_string.split(',')]

    try:
        contrast_col = contrasts[0]
        contrast1 = contrasts[1]
        contrast2 = contrasts[2]

    except Exception as e:
        log.error("error {}".format(e))
        raise SystemExit(e)

    # check that contrast column found in metadata
    assert contrast_col in agg_df.columns

    # check that contrast values are found in the contrast column
    assert contrast1 in agg_df[contrast_col].values
    assert contrast2 in agg_df[contrast_col].values


def compare_two_file_length(file1, file2):
    file_array = "({} {})".format(file1, file2)
    compare_cmd = r'''file_array={}; for file in "${{file_array[@]}}"; do wc -l $file | cut -f1 -d' '; done'''.format \
        (file_array)
    results = subprocess_cmd(compare_cmd, shell=True)
    return results

def check_master_dhs(master_dhs_list):
    if (master_dhs_list.lower() != 'none') and (master_dhs_list is not None):
        if os.path.exists(master_dhs_list) and os.stat(master_dhs_list).st_size != 0:

            log.info("Master DHS file: {} ".format(master_dhs_list))
        else:
            log.error("Master DHS file {} does not exist or is empty.".format(master_dhs_list))
            raise SystemExit("Master DHS file {} does not exist or is empty.".format(master_dhs_list))
    else:
        log.info("Master DHS file will be created.")

def check_agg_list(agg_df):

    # check that input agg_list contains agg_id's and a condition column
    if (len(agg_df['agg_id']) < 1) or (len(agg_df['condition']) < 1):
        log.error("{} must contain columns titled 'agg_id' and 'condition'".format(agg_df))
        raise SystemExit("{} must contain columns titled 'agg_id' and 'condition'".format(agg_df))

def check_exp_design(experimental_design, input_contrasts, agg_df):

    # create list of terms in experimental design
    exp_design = [x.replace('"', '').strip(" ") for x in experimental_design.replace("~", "").split("+")]

    for col in exp_design:
        # make sure design is not empty
        if exp_design is not None:

            # make sure that exp.design term is a column in input metadata
            if col in agg_df.columns:
                log.info("{} column from experimental design found in input metadata.".format(col))
                pass
            else:
                print(col)
                log.error("{} column not found in input metadata".format(col))
                raise SystemExit("Column '{}' not found in input metadata".format(col))
        else:
            log.error("Experimental design cannot be left blank or equal to none.")
            raise SystemExit("Experimental design cannot be left blank or equal to none.")

    # parse experimental design variable
    if ((len(exp_design) < 1) or (exp_design == "none")):
        log.error("Experimental design must be entered. Ex: 'taxonomy' or 'taxonomy + treatment'")
        raise SystemExit("Experimental design must be entered. Ex: 'taxonomy' or 'taxonomy + treatment'")
    # if experimental design contains more than one variable, contrasts are required
    elif len(exp_design) > 1:
        print(exp_design)
        if len(input_contrasts.split(',')) < 3:
            log.error("When using more than one experimental design parameter, you must specify contrasts.")
            raise SystemExit("When using more than one experimental design parameter, you must specify contrasts.")
    else:
        log.info("Experimental design: ~{}".format(experimental_design))

def check_wd(working_dir):
    if os.path.exists(working_dir):
        log.info("Working directory set to: {}".format(working_dir))
    else:
        log.info("Working directory {} does not exist. Creating.".format(working_dir))
        os.makedirs(working_dir)

def check_sample_threshold(sample_threshold, agg_df):
    agg_count = len(agg_df)
    if (int(sample_threshold) == 0):
        log.info("Samples will not be filtered".format(sample_threshold))
    elif (int(sample_threshold) > 0):
        if int(sample_threshold) < agg_count:
            log.info("Regions that do not occur in more than {} samples will be filtered out.".format(sample_threshold))
        else:
            log.error("Sample threshold must be an integer less than the number of samples entered.")
            raise SystemExit("Sample threshold must be an integer less than the number of samples entered.")
    else:
        log.error("Sample threshold must be entered as an integer. Set to 0 for no filtering.")
        raise SystemExit("Sample threshold must be entered as an integer. Set to 0 for no filtering.")

def check_norm_method(norm_method):
    if norm_method.lower() == 'custom':
        log.info("Custom normalization factors created using John Lazar's normalization script.")
    elif (norm_method.lower() == 'none') or (norm_method is None):
        log.info("Normalization factors will be automatically calcuated in R.")
    else:
        log.error("Normalization method must be either 'custom' or 'none'.")
        raise SystemExit("Normalization method must be either custom or none.")

def check_geomean(geomean):
    if geomean.lower() == 'true':
        log.info("Normalization method = geometric mean.")
    elif geomean.lower() == 'false':
        log.info("Normalizatoin method = Loess group norm")
    else:
        log.error("Please specify true to use geomean or false to use Loess group norm.")
        raise SystemExit("Please specify true to use geomean or false to use Loess group norm.")

def check_combined_counts(combined_counts_file, agg_list):
    count_combined_cols_cmd = "head -n 1 {} | awk '{{print NF}}'".format(combined_counts_file)
    col_count = subprocess_cmd(count_combined_cols_cmd, shell=True)
    agg_count_cmd = "cat {} | wc -l".format(agg_list)
    agg_count = subprocess_cmd(agg_count_cmd, shell=True)
    if int(col_count[0]) != int(agg_count[0]):
        log.error("Number of colums in combined counts file does not match input aggregate id's.")
        raise SystemExit("Number of colums in combined counts file does not match input aggregate id's.")

def check_label_col(label_col, agg_df):
    if label_col.lower() != "condition":
        try:
            assert label_col in agg_df.columns
        except AssertionError:
            log.error("Label column not found in metadata file")
            raise SystemExit("Label column not found in metadata file")

def check_metadata_nulls(agg_df):
    # raise SystemExit if input metadata contains nulls
    if agg_df['agg_id'].isnull().values.any():
        log.error(agg_df.isnull().sum())
        raise  SystemExit("Input metadata contains nulls in the following column(s): \n {}".format(agg_df.isnull().sum()))

def check_duplicates(agg_df):
    # check that number of unique agg_id's is equal to number of input agg_id's (no dups)
    assert len(agg_df.agg_id.unique()) == len(agg_df['agg_id'])

def check_symlinks(agg_df, wd):
    symlink_dir = "{}/input_files/".format(wd)

    # count number of symlinks to hotspots and cutcounts files
    cutcounts = []
    base = []

    # list input directory contents
    for file in os.listdir(symlink_dir):
        # create list of cutcounts files
        if file.endswith("cutcounts.bed.starch"):
            cutcounts.append(file)
        # create list of hotspots files
        if file.endswith("hotspots.bed.starch") or file.endswith("peaks.bed.starch"):
            base.append(file)

    # count up hotspots, cutcounts and input metatdata (minus header)
    cutcount_file_count = len(cutcounts)
    hotspots_file_count = len(base)
    # calculate number of agg_id's in metadata minus header
    metadata_sample_count = (len(agg_df))
    print("{} cutcount files symlinked.".format(cutcount_file_count))

    print("{} hotspot or peaks files symlinked.".format(hotspots_file_count))
    print("{} samples in metadata file.".format(metadata_sample_count))
    # check that there are equal counts and hotspots files
    assert cutcount_file_count == hotspots_file_count

    # check that the number of files is same length as agg_df - header
    assert metadata_sample_count == hotspots_file_count

def main(master_dhs_list, experimental_design,
         working_dir, sample_threshold, norm_method, geomean,
         contrasts, label_col, agg_df):

    check_wd(working_dir=working_dir)
    check_agg_list(agg_df=agg_df)
    check_master_dhs(master_dhs_list= master_dhs_list)
    check_contrasts(input_contrast_string=contrasts, agg_df=agg_df)
    check_exp_design(experimental_design= experimental_design, input_contrasts=contrasts, agg_df=agg_df)
    check_sample_threshold(sample_threshold= sample_threshold, agg_df=agg_df)
    check_norm_method(norm_method= norm_method)
    check_geomean(geomean= geomean)
    check_label_col(label_col=label_col, agg_df=agg_df)
    check_metadata_nulls(agg_df=agg_df)
    check_duplicates(agg_df=agg_df)

if __name__ == "__main__":
    main()


