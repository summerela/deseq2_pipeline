import fileinput
import os
import pdb
import shutil
import subprocess
import sys
from collections import OrderedDict

import gridmap
import numpy as np

from deseq_pipeline.deseq.normalization import bedtools

BEDOPS = 'bedops'#'/home/jlazar/.local/bedops/bedops'
BEDMAP = 'bedmap'#'/home/jlazar/.local/bedops/bedmap'
UNSTARCH = '/home/jlazar/.local/bedops/unstarch'
MOTIF_MAPPINGS = '/home/jlazar/scripts/final.gene2motif-mappings.v2.hg19.bed'
CORRELATION_FILE = '/net/lebowski/vol1/work/erynes/data/GenomicsOfGeneRegulation/correlations_37celltypes_inclKBM7K562andSCLC_alsoKTandRC/corrs_distalsFirst_above0.7_37celltypes.bed8'

class Feature_List:
    def __init__(self,file_list,file_location,
                bed_file = None,count_file = None):
        self.file = file_location
        self.exists = os.path.exists(self.file)
        self.base_folder = os.path.split(self.file)[0]  
        self.typeID = 0    
        if self.exists and self.file:
            self.from_file()
        elif bed_file != None:
            lines = open(bed_file).readlines()
            names = [l.split()[3] for l in lines]
            self.features = OrderedDict()
            count = 0
            for n in names:
                if n not in self.features:
                    self.features[n] = count
                    count += 1
            self.names = [n for n in self.features]
            self.save_list()
        elif count_file != None:
            lines = open(count_file).readlines()
            names = [l.split()[0] for l in lines]
            self.features = OrderedDict()
            count = 0
            for n in names:
                if n not in self.features:
                    self.features[n] = count
                    count += 1
            self.names = [n for n in self.features]
            self.save_list()        
        else:
            self.create_feature_list(file_list)
            self.save_list()
        self.element_number = len(self.features.keys())
    
    def write_subset(self,subset,outfile,genes = False):
        count = 0
        wrote = 0
        subset = set(subset)
        f2 = open(outfile,'w')
        with open(self.file) as f1:
            for line in f1:
                if count in subset:
                    if not genes:
                        f2.write(line.strip()[0]+'\n')
                    else:
                        eng = line.strip().split()[0]
                        line = genes.dict[eng].attr['gene_name'] + '\n'
                        f2.write(line)
                        wrote += 1
                count += 1
        f2.close()       
    
    def from_file(self):
        self.features = OrderedDict({})
        self.names = []
        with open(self.file) as f1:
            for line in f1:
                f,index = line.strip().split('\t')
                self.features[f] = int(index)
                self.names.append(f)
    
    def create_feature_list(self,file_list):
        count = 0
        features = OrderedDict({})
        self.names = []
        for file in file_list:
            header = True
            f1 = open(file)
            for line in f1:
                if header:
                    header = False
                    continue
                data = line.split()
                if data[0] not in features:
                    features[data[0]] = count
                    count += 1
            f1.close()    
        self.features = features
        self.names = [k for k in features]
    def save_list(self):
        with open(self.file,'w') as f1:
            f1.write('\n'.join([f+'\t'+str(self.features[f]) for f in self.features]))
    def associate_with_annotation(self,annotation,feature_type = 'gene_id',get_discrepancy = False): 
        ### right now this is essentially only ensembID to gene name conversion
        self.gene_association = {}
        ref_not_in_gene = []
        gene_not_in_ref = []
        annotation.make_attribute_dict(feature_type)
        annotation_dict = getattr(annotation,feature_type+'_dict')
        for r in self.features:
            if r in annotation.dict:
                self.gene_association[annotation.dict[r].attr[feature_type]] = self.features[r]
            else: 
                new = r[0:r.find('.')]
                if new in annotation_dict:
                    self.gene_association[annotation_dict[new]] = self.features[r]
                elif get_discrepancy:
                    ref_not_in_gene.append(new)
        if get_discrepancy:
            gene_not_in_ref = [g for g in annotation_dict if g not in self.gene_association]
        return ref_not_in_gene,gene_not_in_ref
        ##### need to incorporate isoforms into this/better genbank parsing            
    
class DHS_List:
    def __init__(self,file_location):
        self.file = file_location
        self.exists = os.path.exists(self.file)
        self.base_folder = os.path.split(self.file)[0]  
        self.typeID = 0  
        
    def get_locations(self):
        import HTSeq
        bed_reader = HTSeq.BED_Reader(self.file)
        self.locations = [r for r in bed_reader]
        
    def write_subset(self,subset,outfile,values = [],genes = None):
        count = 0
        wrote = 0
        subset = set(subset)
        f2 = open(outfile,'w')
        with open(self.file) as f1:
            for line in f1:
                if count in subset:
                    if len(values) != len(subset):
                        f2.write(line)
                    else:
                        line = line.strip() + '\t' + str(values[wrote])+'\n'
                        f2.write(line)
                        wrote += 1
                count += 1
        f2.close()        
    
    def correlated_tss(self,tss_file = CORRELATION_FILE,correlation_limit = 0.7):  ### correlations limit is set to eric's default
        cmd = subprocess.Popen(['bedmap','--echo','--echo-map-id',self.file,tss_file],stdout=subprocess.PIPE)
        self.corr_gene_ids = ['']*self.element_number
        count = 0
        for line in iter(cmd.stdout.readline,''): 
            id,names = line.split('|')
            self.corr_gene_ids[count] = list(set(names.strip().split(';')))
            count += 1
        cmd.stdout.close()
    
    def closest_tss(self,gene_annotation,promoter_limit = 1000,distal_limit = 250000,tolerance = 'strict'):
        if tolerance == 'strict':
            tss_file = gene_annotation.RefSeq_tss_file
        else:
            tss_file = gene_annotation.tss_file
        cmd = subprocess.Popen(['closest-features','--closest','--dist','--no-ref',self.file,tss_file],stdout=subprocess.PIPE)
        self.gene_ids = ['']*self.element_number
        self.distances = np.zeros(self.element_number)
        ##### {intergenic = 0, promoter = 1, distal = -1}
        self.promoter = np.zeros(self.element_number,dtype = bool)
        self.distal = np.zeros(self.element_number,dtype = bool)
        count = 0
        for line in iter(cmd.stdout.readline,''): 
            gene,distance = line.split('|')
            self.distances[count] = int(distance)
            gene_name = gene.split()[3]
            self.gene_ids[count] = gene_name
            count += 1
        cmd.stdout.close()
        promoter = np.where(np.abs(self.distances) < promoter_limit)[0]
        self.promoter[promoter] = True
        within_range = np.where(np.abs(self.distances) < distal_limit)[0]
        distal = np.setdiff1d(within_range,promoter)
        self.distal[distal] = True
    def tss_within_distance(self,gene_annotation,distance = 500000):
        f1 = open('/home/jlazar/tmp/error.txt','w')
        cmd = subprocess.Popen(['bedmap','--echo','--echo-map-id','--range',
                                str(distance),self.file,gene_annotation.tss_file],
                                stdout = subprocess.PIPE,stderr = f1)
        self.associated_tss = []
        for line in iter(cmd.stdout.readline,''):
            genes = line.strip().split('|')[1].split(';')
            self.associated_tss.append(list(set(genes)))
        cmd.stdout.close()
        f1.close()
        
    
    def create_motif_file(self):
        self.motif_folder = os.path.join(self.base_folder,'motifs')
        self.motif_file = os.path.join(self.motif_folder,'all.motifs'+str(self.typeID))
        if os.path.exists(self.motif_file):
            return
        unstarch_x = [UNSTARCH,'/home/jvierstra/data/motifs/fimo/hg19.human.xfac.1e-4/fimo.combined.xfac.1e-4.parsed.starch']
        unstarch_j = [UNSTARCH,'/home/jvierstra/data/motifs/fimo/hg19.human.jaspar.1e-4/fimo.combined.jaspar.1e-4.parsed.starch']
        unstarch_u = [UNSTARCH,'/home/jvierstra/data/motifs/fimo/hg19.human.uniprobe.1e-4/fimo.combined.uniprobe.1e-4.parsed.starch']        
        combine = [BEDOPS,'-e','-1','-',self.file]
        if not os.path.exists(self.motif_folder):
            os.mkdir(self.motif_folder)
        xfac_location = os.path.join(self.motif_folder,'xfac.motifs'+str(self.typeID))
        jaspar_location = os.path.join(self.motif_folder,'jaspar.motifs'+str(self.typeID))
        uniprobe_location = os.path.join(self.motif_folder,'uniprobe.motifs'+str(self.typeID))
        
        p1 = subprocess.Popen(unstarch_x,stdout = subprocess.PIPE)
        p2 = subprocess.Popen(combine,stdin=p1.stdout,stdout=open(xfac_location,'w'))
        p1.stdout.close()
        p2.communicate()

        p1 = subprocess.Popen(unstarch_j,stdout = subprocess.PIPE)
        p2 = subprocess.Popen(combine,stdin=p1.stdout,stdout=open(jaspar_location,'w'))
        p1.stdout.close()
        p2.communicate()

        p1 = subprocess.Popen(unstarch_u,stdout = subprocess.PIPE)
        p2 = subprocess.Popen(combine,stdin=p1.stdout,stdout=open(uniprobe_location,'w'))
        p1.stdout.close()
        p2.communicate()
 
        
        all_locations = os.path.join(self.motif_folder,'all.motifs'+str(self.typeID))
        all = [BEDOPS,'-u',xfac_location,jaspar_location,uniprobe_location]
        run = subprocess.Popen(all,stdout=open(all_locations,'w'))
        run.wait()         

class Centroid_List(DHS_List):
    def __init__(self,experiment_list,file_location,feature = ''):
        DHS_List.__init__(self,file_location)  
        self.vector_storage = os.path.join(self.base_folder,'master_vectors')
        if not feature:
            self.feature_type = experiment_list[0].feature_data
        else:
            self.feature_type = feature               
        if not os.path.exists(file_location):
            self.create_master_list(experiment_list)
        
        self.overlap_criterion = '--fraction-map'
        self.overlap_number = '0.25'
        self.typeID = 1 
        count = 0
        with open(self.file) as f1:
            for line in f1:
                count += 1
        self.element_number = count
        
    def create_master_list(self,experiment_list):
        peak_files = [exp.data_files[self.feature_type] for exp in experiment_list]
        density_files = [exp.data_files[exp.magnitude_data] for exp in experiment_list]
        #merge_job = gridmap.Job(bedtools.bedops,[peak_files,'tempfile',['-m']])####
        #merge_job.execute()
        merge_file = bedtools.bedops(peak_files, 'tempfile', ['-m'])
        fold = open(merge_file)
        with open(self.file,'w') as f1:
            count = 0
            for line in fold:
                f1.write(line.strip()+'\t%d\n' %count)
                count += 1
        fold.close()        
        os.remove(merge_file)
        return
        shutil.move(merge_file,self.file)
        return ### start out with just a merge file
        partition_job = gridmap.Job(bedtools.bedops([merge_file], 'tempfile', ['--chop', '20']))####
        partition_job.execute()        
        partition_file = partition_job.ret
        mapping_jobs = [gridmap.Job(bedtools.bedmap([partition_file, d], 'tempfile', ['--mean', '--exact', '--faster']))
                        for d in density_files]
        mapping_files = gridmap.process_jobs(mapping_jobs)
        with open(partition_file) as partition:
            previous = ''
            mapping_file_objects = [open(f) for f in mapping_files]
            merge_file_object = open(merge_file)
            final_file_object = open(self.file,'w')
            count = 0
            merge_element = merge_file_object.readline()
            for line in partition:
                data = line.split.strip()
                values = [float(r.readline().strip()) for r in mapping_file_objects]
                if data[1] == previous:   # this should be based on merge_element
                    scores.append(sum(values))
                else:
                    new_elements = parse_peak_scores(scores)
                    if not new_elements:
                        final_file_object.write(merge_element)
                    #### need to write the elements here
                    scores = [sum(values)]
                    count += 1
                    merge_element = merge_file_object.readline()
                previous = data[2]
        #### then we should partition into 20bp windows to look at the distribution
        #### could do 99% of the total mass to narrow
        #### test each of bimodality?
def parse_peak_scores(density_list,trim_threshold=.995):
    lmaxs = [density_list[i]<density_list[i+1]>density_list[i+2] for i in xrange(len(density_list-2))]
    total_density = sum(density_list)
    if sum(lmaxs) < 2:
        kept = trim_density(densities,total_density,trim_threshold)
    return kept
def trim_density(density_list,total_density,trthreshold):
    max_trim = total_density*(1-trim_threshold)
    edges = [0,-1]
    trimmed = min(density_list[left[0]],density_list[left[-1]])
    while trimmed < max_trim:
        if density_list[left[0]]<density_list[left[-1]]:
            edges[0] += 1
        else:
            edges[1] += 1
        trimmed += min(density_list[left[0]],density_list[left[-1]])
    
    
        
        
class Master_List(DHS_List):

    def __init__(self,experiment_list,file_location,feature = ''):
        DHS_List.__init__(self,file_location)
        self.vector_folder = os.path.join(self.base_folder,'master_vectors')
        if not os.path.exists(self.vector_folder):
            os.mkdir(self.vector_folder)
        if not feature:
            self.feature_type = experiment_list[0].feature_data
        else:
            self.feature_type = feature               
        self.create_master_list(experiment_list)
        self.overlap_criterion = '--fraction-map'
        self.overlap_number = '0.25'
        self.typeID = 1
    
    def create_master_list(self,experiment_list,overwrite = False):
        if not overwrite and self.exists:
            self.element_number = 0
            with open(self.file) as f1:
                for line in f1: self.element_number += 1
            return
        else:
            fob = open(self.file,'w')
            fob.close()
        self.union_file = os.path.join(self.base_folder,'union_file.bed')
        if overwrite or not os.path.exists(self.union_file):
            union_temp = self.union_file + '.temp'
            f1 = open(self.union_file,'w')
            f1.close()   
            for exp in experiment_list:
                os.rename(self.union_file,union_temp)
                with open(self.union_file,'w') as fOb:
                    #extract = ['gunzip','-c',exp.data_files[self.feature_type]]
                    extract = ['gunzip','-cf',exp.data_files[self.feature_type]]
                    gunzip = subprocess.Popen(extract,stdout=subprocess.PIPE)
                    cut = ['cut','-f1-4,8']
                    cutting = subprocess.Popen(cut,stdout=subprocess.PIPE,stdin=gunzip.stdout)
                    cmd = [BEDOPS ,'-u',union_temp,'-']
                    union = subprocess.Popen(cmd,stdout=subprocess.PIPE,stdin=cutting.stdout)
                    uniq_cmd = ['uniq']
                    uniq = subprocess.Popen(uniq_cmd,stdout=fOb,stdin=union.stdout)
                    gunzip.stdout.close()
                    union.stdout.close()
                    uniq.communicate()
                os.remove(union_temp)   
        self.union_file_subset = os.path.join(self.base_folder,'remaining_dhs.bed')
        shutil.copyfile(self.union_file,self.union_file_subset)
        fOb = open(self.file + '.mid','w')
        fOb.close()
        while os.stat(self.union_file_subset).st_size !=0:
            self.master_list_subroutine(self)
        os.remove(self.union_file_subset)
        os.remove(self.union_file_subset+'.temp')
        os.remove(self.file + '.temp')
        with open(self.file,'w') as f1:
            cutting = subprocess.Popen(['cut','-f1-3',self.file + '.mid'],stdout = subprocess.PIPE)
            uniquing = subprocess.Popen(['uniq'],stdin = cutting.stdout,stdout = f1)
            cutting.stdout.close()
            uniquing.communicate()
        os.remove(self.file+'.mid')
        count = 0
        for line in fileinput.input(self.file,inplace=True):
            line = line.strip()+'\t'+str(count)+'\n'
            sys.stdout.write(line)
            count += 1
        self.element_number = count
        os.remove(self.union_file)
        
    def master_list_subroutine(self,original):    
        os.rename(self.file+'.mid',self.file+'.temp')
        merge_cmd = [BEDOPS,'-m',self.union_file_subset]
        merge = subprocess.Popen(merge_cmd,stdout=subprocess.PIPE)
        map_cmd = [BEDMAP,'--max-element','-',self.union_file_subset]
        map = subprocess.Popen(map_cmd,stdout=subprocess.PIPE,stdin=merge.stdout)
        union_cmd = [BEDOPS,'-u','-',self.file+'.temp']
        with open(self.file+'.mid','w') as f1:
            union = subprocess.Popen(union_cmd,stdout=f1,stdin=map.stdout)
            merge.stdout.close()
            map.stdout.close()
            union.communicate()
        os.rename(self.union_file_subset,self.union_file_subset+'.temp')
        leftover_cmd = [BEDOPS,'-n','-25%',self.union_file_subset+'.temp',self.file+'.mid']
        with open(self.union_file_subset,'w') as f1:
            leftover = subprocess.Popen(leftover_cmd,stdout=f1)
            leftover.communicate()    

    def get_overlaps(self,second_file):
        ovr = bedtools.bedmap([self.file, second_file], 'pipe', ['--indicator',
                                                                 self.overlap_criteria, self.overlap_number])
        return np.loadtxt(ovr)
            
class Merge_List(DHS_List):

    def __init__(self,experiment_list,file_location):
        DHS_List.__init__(self,file_location)
        self.vector_storage = os.path.join(self.base_folder,'union_vectors')
        self.create_union_vector(experiment_list)
        self.overlap_criterion = '--fraction-map'
        self.overlap_number = '1'
        self.typeID = 2

    def create_union_vector(self,experiment_list, overwrite = False):
        if not overwrite and self.exists:
            self.element_number = 0
            with open(self.file) as f1:
                for line in f1: self.element_number += 1
            return
        else:
            fob = open(self.file,'w')
            fob.close()
        self.union_file = os.path.join(self.base_folder,'union_file.bed')
        union_temp = self.union_file + '.temp'
        f1 = open(self.union_file,'w')
        f1.close()
        for exp in experiment_list:
            os.rename(self.union_file,union_temp)
            with open(self.union_file,'w') as fOb:
                extract = ['gunzip','-c',exp.data_files['peaks']]
                gunzip = subprocess.Popen(extract,stdout=subprocess.PIPE)
                cmd = [BEDOPS ,'-u',union_temp,'-']
                union = subprocess.Popen(cmd,stdout=subprocess.PIPE,stdin=gunzip.stdout)
                uniq_cmd = ['uniq']
                uniq = subprocess.Popen(uniq_cmd,stdout=fOb,stdin=union.stdout)
                gunzip.stdout.close()
                union.stdout.close()
                uniq.communicate()
            os.remove(union_temp)
            
        with open(self.file,'w') as f2:
            cmd = [BEDMAP, '--echo-map-range', '--fraction-map', '0.25', self.union_file, self.union_file]
            ref = subprocess.Popen(cmd,stdout=subprocess.PIPE)
            uniq = subprocess.Popen(['uniq'],stdout=f2,stdin=ref.stdout)             
            ref.stdout.close()
            uniq.communicate()
        
        count = 0            
        for line in fileinput.input(self.file,inplace=True):
            line = line.strip()+'\t'+str(count)+'\n'
            sys.stdout.write(line)
            count += 1
        self.element_number = count            
        
        
        
class Master_List_Outlier(DHS_List):

    def __init__(self,experiment_list,file_location):
        DHS_List.__init__(self,file_location)
        self.vector_storage = os.path.join(self.base_folder,'master_vectors')
        self.create_master_list(experiment_list)
        self.overlap_criterion = '--fraction-map'
        self.overlap_number = '0.25'
        self.typeID = 1
    
    def create_master_list(self,experiment_list,overwrite = False):
        if not overwrite and self.exists:
            self.element_number = 0
            with open(self.file) as f1:
                for line in f1: self.element_number += 1
            return
        else:
            fob = open(self.file,'w')
            fob.close()
        self.union_file = os.path.join(self.base_folder,'union_file.bed')
        if not overwrite or not os.path.exists(self.union_file):
            union_temp = self.union_file + '.temp'
            f1 = open(self.union_file,'w')
            f1.close()   
            for exp in experiment_list:
                os.rename(self.union_file,union_temp)
                with open(self.union_file,'w') as fOb:
                    extract = ['gunzip','-c',exp.data_files['peaks']]
                    gunzip = subprocess.Popen(extract,stdout=subprocess.PIPE)
                    cut = ['cut','-f1-4,8']
                    cutting = subprocess.Popen(cut,stdout=subprocess.PIPE,stdin=gunzip.stdout)
                    cmd = [BEDOPS ,'-u',union_temp,'-']
                    union = subprocess.Popen(cmd,stdout=subprocess.PIPE,stdin=cutting.stdout)
                    uniq_cmd = ['uniq']
                    uniq = subprocess.Popen(uniq_cmd,stdout=fOb,stdin=union.stdout)
                    gunzip.stdout.close()
                    cmd.stdout.close()
                    union.stdout.close()
                    uniq.communicate()
                os.remove(union_temp)   
        self.union_file_subset = os.path.join(self.base_folder,'remaining_dhs.bed')
        shutil.copyfile(self.union_file,self.union_file_subset)
        while os.stat(self.union_file_subset).st_size !=0:
            self.master_list_subroutine(self)
        os.remove(self.union_file_subset)
        os.remove(self.union_file_subset+'.temp')
        os.remove(self.file + '.temp')
        count = 0    
        for line in fileinput.input(self.file,inplace=True):
            line = line.strip()+'\t'+str(count)+'\n'
            sys.stdout.write(line)
            count += 1
        self.element_number = count
        
    def master_list_subroutine(self,original):    
        os.rename(self.file,self.file+'.temp')
        merge_cmd = [BEDOPS,'-m',self.union_file_subset]
        merge = subprocess.Popen(merge_cmd,stdout=subprocess.PIPE)
        map_cmd = [BEDMAP,'--echo-map','-',self.union_file_subset]
        map = subprocess.Popen(merge_cmd,stdout=subprocess.PIPE,stdin=merge.stdout)
        union_cmd = [BEDOPS,'-u','-',self.file+'.temp']
        with open(self.file,'w') as f1:
            union = subprocess.Popen(union_cmd,stdout=f1,stdin=map.stdout)
            merge.stdout.close()
            map.stdout.close()
            union.communicate()
        pdb.set_trace()
        os.rename(self.union_file_subset,self.union_file_subset+'.temp')
        leftover_cmd = [BEDOPS,'-n','-25%',self.union_file_subset+'.temp',self.file]
        with open(self.union_file_subset,'w') as f1:
            leftover = subprocess.Popen(leftover_cmd,stdout=f1)
            leftover.communicate()          
