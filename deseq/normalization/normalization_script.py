#! /usr/bin/env python
import argparse
import logging
import os
import time
import numpy as np
import normalization

# script log created as deseq2_{today}.log in working directory
today = time.strftime("%d%m%Y")
FORMAT = "[%(name)s | Line:%(lineno)s| Function:%(funcName)s ] %(message)s"
logging.basicConfig(filename="deseq2_{}.log".format(today), level=logging.INFO, filemode='w',
                    format=FORMAT)
log = logging.getLogger('normalization_script')

##### simple script to run the normalization
#### takes a file with a table of counts as the first argument
### and a table of peak calls as the second
### outputs the size factors to use with DESeq


def run_normalization_file(count_file,calls_file,geomean,output_file = None):
    density_data = np.loadtxt(count_file)
    calls_data = np.loadtxt(calls_file)
    calls_data = calls_data.astype(bool)
    #print density_data.shape
    norm_number = [1/(np.sum(density_data[:,i])*2/100000.) for i in range(density_data.shape[1])]
    
    density_vectors = {i: density_data[:,i]*norm_number[i] for i in range(density_data.shape[1])}
    calls_vectors = {i: calls_data[:,i] for i in range(calls_data.shape[1])}

    if geomean:
        normed_vectors = normalization.loess_geomean_norm(density_vectors, calls_vectors)
    else:
        normed_vectors = normalization.loess_group_norm(density_vectors, calls_vectors)
    base_name = os.path.split(count_file)[1].strip('.txt')
    base_dir = os.path.split(count_file)[0]

    orders = range(calls_data.shape[1])
    if output_file == None:
        output_file_factors = os.path.join(base_dir,base_name + '.size_factors')
        output_file_counts = os.path.join(base_dir,base_name + '.norm_counts')
    else:
        output_file_factors = output_file.rstrip('.txt') + '.size_factors'
        output_file_counts = output_file.rstrip('.txt') + '.norm_counts'
    
    density_vectors = {i:density_data[:,i] for i in orders}
    normed_vectors = {i:normed_vectors[i]/np.mean(norm_number) for i in orders}
    normalization_file = normalization.create_normalization_file(normed_vectors,
                                                                 density_vectors,
                                                                 orders, normalization_file = output_file_factors) #### need to get a non-random file
    np.savetxt(output_file_counts,np.transpose(np.array([normed_vectors[r] for r in orders])),fmt='%.6e')  

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process Count Files for Normalization.')
    parser.add_argument('counts_file',help = 'file of counts')
    parser.add_argument('calls_file',help = 'file of peak calls')
    parser.add_argument('--geomean',help='normalize to geometric mean',
                        default = False,action='store_true')
    parser.add_argument('-o','--output',type = str,help='output file',
                        default = None)                        
    
    
    args = parser.parse_args()

    try:
        run_normalization_file(args.counts_file,args.calls_file,args.geomean,args.output)

    except Exception as e:
        log.error(e)
