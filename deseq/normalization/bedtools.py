import numpy as np
import os, subprocess, multiprocessing
import random,string

try:
    import gridmap
    import HTSeq
except: pass

TAB = '\t'
ENDLINE = '\n'
JARCH_LOCATION = '/home/rthurman/dev/CompBio/utility/bin/gchr'
#TEMP_FOLDER = '/home/jlazar/tmp'
TEMP_FOLDER = os.getcwd()

def return_scores(bed_file):
    return np.loadtxt(bed_file,usecols=(4,))

def combine_bed_files(bed_list,out_file='tmpfile'):
    all = []
    count = 1
    intermed = make_tmpfile()
    if out_file == 'tmpfile':
        out_file = make_tmpfile()
    for f in bed_list:
        f2 = open(f)
        for line in f2:
            all.append('\t'.join(line.strip().split(TAB)[0:3])+TAB+'.'+TAB+str(count))
        count += 1
        f2.close()
    with open(intermed,'w') as f1:  
        f1.write('\n'.join(all))
    sort = ['sort-bed',intermed]
    with open(out_file,'w') as f2:
        run = subprocess.Popen(sort,stdout=f2)
        run.communicate()
    os.remove(intermed)   
    return out_file
    
def parse_bed(file):
    elements = []
    fOb = open(file)
    for line in fOb:
        data = line.strip().split('\t')
        elements.append(structure(data[0:6]))
    fOb.close()   
    return elements
    
def count_lines(file):
    count = 0
    with open(file) as f1:
        for line in f1:
            if line:
                count += 1
    return count

def make_bed(feature_list,file='tmpfile',structures = False):
    if structures:
        features = [c.chr,c.start,c.end,c.name]
    else:
        features = feature_list
    lines = ['\t'.join(f) for f in features]
    if file == 'tmpfile':
        file = make_tmpfile()
    with open(file,'w') as f1:
        f1.write('\n'.join(lines))
    return file
    
class structure:
    def init(self,chr,start,end,name,strand):
        self.chr = chr
        self.start = start
        self.end = end
        self.name = name
        self.strand = strand
      
def on_cluster(functions,function_args,tmp_dir='/home/jlazar/tmp'):
    all_jobs = []
    for i,arg in enumerate(function_args):
        cluster_job = gridmap.Job(functions[i],arg)
        all_jobs.append(cluster_job)
    gridmap.process_jobs(all_jobs,temp_dir=tmp_dir)

def more_than_one(input_files,output_file):    
    if output_file == 'tempfile':
        output_file = make_tmpfile()
    parsed_files = []
    for file in input_files:
        parsed_files.append(parse_input_file(file))        
    tmpf = make_tmpfile()
    with open(tmpf,'w') as f1:
        merge = subprocess.Popen(['bedops','-u']+[p[0] for p in parsed_files],stdout=f1)
        merge.communicate()
    rsym = subprocess.Popen(['bedmap','--count',tmpf,tmpf],stdout=subprocess.PIPE)            
    count = 0
    lines = []          
    f_base = open(tmpf)
    for line in iter(rsym.stdout.readline,''):
        print_line = f_base.readline()
        if int(line.strip()) > 1:
            lines.append(print_line)
    f_base.close()
    rsym.stdout.close()    
    with open(output_file,'w') as f2:
        f2.write(''.join(lines))
    for filename,tempfile in parsed_files:
        if tempfile:
            os.remove(filename)
    os.remove(tmpf)
    return output_file
    
def bedops(input_files,output_file,args):
    if output_file == 'tempfile':
        output_file = make_tmpfile()
    parsed_files = []
    for file in input_files:
        parsed_files.append(parse_input_file(file))        
    with open(output_file,'w') as f1:
        if output_file != 'pipe':
            rbedmap = subprocess.Popen(['bedops']+args+[p[0] for p in parsed_files],stdout=f1)
            rbedmap.communicate()
        else:
            rbedmap = subprocess.Popen(['bedops']+args+[p[0] for p in parsed_files],stdout=subprocess.PIPE)
            output_file = rbedmap.stdout

    for filename,tempfile in parsed_files:
        if tempfile:
            os.remove(filename)
    return output_file

def distance_beds(associations_dict,bed1,bed2):
    ### associations dict is dictionary of names from first with list of secont
    bed1read = HTSeq.BED_Reader(bed1)
    reverse_set = set()
    for k in associations_dict:
        for match in associations_dict[k]:
            reverse_set.add(str(match))
    distances = {} #defaultdict(list)
    second_bed_elements = {}
    bed2read = HTSeq.BED_Reader(bed2)
    for ft2 in bed2read:
        if ft2.name in reverse_set:
            second_bed_elements[ft2.name] = ft2
    for ft1 in bed1read:
        if ft1.name not in associations_dict:
            continue
        distances[ft1.name] = ['' for x in range(len(associations_dict[ft1.name]))]
        for i,match in enumerate(associations_dict[ft1.name]):
            try:
                match_element = second_bed_elements[str(match)]
                if match_element.iv.chrom == ft1.iv.chrom:
                    distance = match_element.iv.start - ft1.iv.start
                    distances[ft1.name][i] = distance
                else:
                    distances[ft1.name][i] = 'trans'            
            except KeyError:
                distances[ft1.name][i] = 'nan'
    return distances
def bedmap(input_files,output_file,args):
    ### might be slightly faster if always pipe the second file
    if output_file == 'tempfile':
        output_file = make_tmpfile()
    parsed_files = []
    for file in input_files:
        parsed_files.append(parse_input_file(file))
    with open(output_file,'w') as f1:
        if output_file != 'pipe':
            rbedmap = subprocess.Popen(['bedmap']+args+[p[0] for p in parsed_files],stdout=f1)
            rbedmap.communicate()
        else:
            rbedmap = subprocess.Popen(['bedmap']+args+[p[0] for p in parsed_files],stdout=subprocess.PIPE)
            output_file = rbedmap.stdout

    for filename,tempfile in parsed_files:
        if tempfile:
            os.remove(filename)
    return output_file
        
def get_file_intensity(feature_file,intensity_file,output_file,score='max'):
        temp = False
        if intensity_file.endswith('.jarch'):
            random_id  = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(15))
            TEMP_FILE = os.path.join(TEMP_FOLDER,random_id+'dhs_tempfile.tmp')            
            temp = True
            density_file = TEMP_FILE
            with open(density_file,'w') as f1:
                unjarch = [JARCH_LOCATION,density_file]
                r_jarch = subprocess.Popen(unjarch,stdout=f1)
                r_jarch.communicate()
        else:
            density_file = intensity_file

        cut = ['/home/jlazar/peak_normalized-density/make_four.sh']
        uncompress = ['gunzip','-cf',feature_file]
        run_uncompress = subprocess.Popen(uncompress,stdout=subprocess.PIPE)
        run_cut = subprocess.Popen(cut,stdin=run_uncompress.stdout,stdout=subprocess.PIPE)
        
        
        get_scores = ['bedmap','--faster','--delim','\t','--echo','--'+score,
                        '-',density_file]
        
        with open(output_file,'w') as f1:
            run_bmap = subprocess.Popen(get_scores,stdin=run_cut.stdout,stdout=f1)
            run_uncompress.stdout.close()
            run_cut.stdout.close()
            run_bmap.communicate()
        if temp:
            os.remove(TEMP_FILE)

def write_subset(input_file,output_file,subset,add_count = False):
    uncompress = ['gunzip','-cf',input_file]
    run_uncompress = subprocess.Popen(uncompress,stdout=subprocess.PIPE)
    f1 = open(output_file,'w')
    count = 0
    subset = set(list(subset))
    for line in iter(run_uncompress.stdout.readline,''):
        if count in subset:
            if add_count:
                f1.write('\t'.join(line.strip().split()[0:3]+[str(count)]) + '\n')
            else:
                f1.write(line)
        count += 1
    f1.close()
    run_uncompress.stdout.close()
            
def uncompress(input_file,output_file,add=False):
    uncompress = ['gunzip','-cf',input_file]   
    if add:
        run = subprocess.Popen(uncompress,stdout=subprocess.PIPE)
        with open(output_file,'w') as f2:
            count = 0
            for line in iter(run.stdout.readline,''):        
                data = line.split('\t')
                data[3] = str(count)
                count += 1
                f2.write('\t'.join(data)+'\n')
        run.stdout.close()   
    else:
        with open(output_file,'w') as f1:
            run = subprocess.Popen(uncompress,stdout=f1)
            run.communicate()
 
def uncompress_pipe(input_file):
    uncompress = ['gunzip','-cf',input_file]   
    run = subprocess.Popen(uncompress,stdout=subprocess.PIPE)
    return run.stdout
    
def get_shared(file1,file2,overlap_criteria):
    random_id = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(15))
    tmp1 = os.path.join(TEMP_FOLDER,random_id+'tmp1.bed')
    tmp2 = os.path.join(TEMP_FOLDER,random_id+'tmp2.bed')
    u1 = multiprocessing.Process(target=uncompress,args=[file1,tmp1])
    u2 = multiprocessing.Process(target=uncompress,args=[file2,tmp2])   
    u1.start()
    u2.start()
    u1.join()
    u2.join()

    non_file_args = ['--indicator']+overlap_criteria
    tmp1o = os.path.join(TEMP_FOLDER,random_id+'tmp1.output.bed') 
    tmp2o = os.path.join(TEMP_FOLDER,random_id+'tmp2.output.bed')     
    u1 = multiprocessing.Process(target=bedmap,args=[[tmp1,tmp2],tmp1o,non_file_args])
    u2 = multiprocessing.Process(target=bedmap,args=[[tmp2,tmp1],tmp2o,non_file_args])
    u1.start()
    u2.start()
    u1.join()
    u2.join()

    out1 = np.loadtxt(tmp1o,usecols=(0,))
    out2 = np.loadtxt(tmp2o,usecols=(0,))
    for file in [tmp1,tmp2,tmp1o,tmp2o]:
        os.remove(file)
    return out1,out2
    
    
def compare2files(file1,file2,overlap_criteria,score=True):
    random_id = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(15))
    tmp1 = os.path.join(TEMP_FOLDER,random_id+'tmp1.bed')
    tmp2 = os.path.join(TEMP_FOLDER,random_id+'tmp2.bed')
    u1 = multiprocessing.Process(target=uncompress,args=[file1,tmp1,True])
    u2 = multiprocessing.Process(target=uncompress,args=[file2,tmp2,True])
    u1.start()
    u2.start()
    u1.join()
    u1.join()
    
    non_file_args = ['bedmap','--echo-map-id']+overlap_criteria
    tmp1o = os.path.join(TEMP_FOLDER,random_id+'tmp1.output.bed') 
    tmp2o = os.path.join(TEMP_FOLDER,random_id+'tmp2.output.bed')     
    u1 = multiprocessing.Process(target=bedmap,args=[[tmp1,tmp2],tmp1o,non_file_args])
    u2 = multiprocessing.Process(target=bedmap,args=[[tmp2,tmp1],tmp2o,non_file_args])
    u1.start()
    u2.start()
    u1.join()
    u2.join()    
    
    over1 = get_bedmap_elements(tmp1o)
    over2 = get_bedmap_elements(tmp2o)
    for file in [tmp1,tmp2,tmp1o,tmp2]:
        os.remove(file)    
    return over1,over2

def parse_input_file(filename,temp_file = True):
    if filename.endswith('.jarch'):
        random_file = make_tmpfile()
        with open(random_file,'w') as f1:
            unjarch = [JARCH_LOCATION,filename]
            r_jarch = subprocess.Popen(unjarch,stdout=f1)
            r_jarch.communicate() 
        return (random_file,True)
    elif filename.endswith('.starch'):
        return (filename,False)
    else:
        random_file = make_tmpfile()
        uncompress(filename,random_file)
        return (random_file,True)
    
def make_tmpfile():
    random_id = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(15))
    return os.path.join(TEMP_FOLDER,random_id+'tmp.bed')

def get_bedmap_elements(file,sep=';'):
    output = []
    with open(file) as f1:
        count = 0
        for line in f1:
            if line=='NAN':
                output.append(())
            else:
                data = line.strip().split(sep)
                output.append(tuple(data))
            count += 1
    return output
        
def create_bed_window(dir_name,window_size):    
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    file = os.path.join(dir_name,'full-genome.'+str(window_size)+'.bed')
    f1 = open(file+'.tmp','w')
    for chr in CHR_SIZES:
        count = 0
        while count+window_size < CHR_SIZES[chr]:
            f1.write(chr+TAB+str(count)+TAB+str(count+window_size)+ENDLINE)
            count += window_size
        f1.write(chr+TAB+str(count)+TAB+str(CHR_SIZES[chr]+1)+ENDLINE)
    f1.close()
    with open(file,'w') as f2:
        r = subprocess.Popen(['sort-bed',file+'.tmp'],stdout=f2)
        r.communicate()
    os.remove(file+'.tmp')

CHR_SIZES = {}
CHR_SIZES['chr1'] = 249250621
CHR_SIZES['chr2'] = 243199373
CHR_SIZES['chr3'] = 198022430
CHR_SIZES['chr4'] = 191154276
CHR_SIZES['chr5'] = 180915260
CHR_SIZES['chr6'] = 171115067
CHR_SIZES['chr7'] = 159138663
CHR_SIZES['chr8'] = 146364022
CHR_SIZES['chr9'] = 141213431
CHR_SIZES['chr10'] = 135534747
CHR_SIZES['chr11'] = 135006516
CHR_SIZES['chr12'] = 133851895
CHR_SIZES['chr13'] = 115169878
CHR_SIZES['chr14'] = 107349540
CHR_SIZES['chr15'] = 102531392
CHR_SIZES['chr16'] = 90354753
CHR_SIZES['chr17'] = 81195210
CHR_SIZES['chr18'] = 78077248
CHR_SIZES['chr19'] = 59128983
CHR_SIZES['chr20'] = 63025520
CHR_SIZES['chr21'] = 48129895
CHR_SIZES['chr22'] = 51304566
CHR_SIZES['chrX'] = 155270560
CHR_SIZES['chrY'] = 59373566
