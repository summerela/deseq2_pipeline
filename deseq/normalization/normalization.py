import logging
import numpy as np
import os
import statsmodels.api as sm
import time

import bedtools

# script log created as deseq2_{today}.log in working directory
today = time.strftime("%d%m%Y")
FORMAT = "[%(name)s | Line:%(lineno)s| Function:%(funcName)s ] %(message)s"
logging.basicConfig(filename="deseq2_{}.log".format(today), level=logging.INFO, filemode='w',
                    format=FORMAT)
log = logging.getLogger('normalization_script')

TMM_RATIOS = {'DNAse':(0.33,0.66),
              'RNA':(0.20,0.8),
              'ChIP-Seq':(.33,.66)}
FILTER_CALLED = {'DNAse':True,
                 'RNA':True,
                 'ChIP-Seq':True}
PEAK_THRESHOLDS = {'DNAse':(0.05,.995),
               'RNA':(0.05,.99),
              'ChIP-Seq':(0.05,.99)}
TOP_THRESHOLDS = {'DNAse':0.025,
               'RNA':0.025,
              'ChIP-Seq':0.025}    

SEED_NUMBER = 1832245 #np.random.RandomState(189342245)
SEED_ITERATION = 10000
SEED = np.random.RandomState(SEED_NUMBER)
        
def update_seed():
    global SEED_NUMBER
    global SEED    
    
    direction = 1
    if SEED_NUMBER > 4294000000:
        direction = -1
        SEED_ITERATION 
    elif SEED_NUMBER < 100000:
        direction = 1
    SEED_NUMBER += direction*(SEED_ITERATION+SEED.randint(1,SEED_ITERATION))
    SEED.seed(SEED_NUMBER)
              
def tmm_graph(density_vectors,called_vectors,data_type = 'DNAse'):
    normed_values = norm_by_graph(density_vectors,called_vectors,method=tmm_norm_pair,
                    data_type=data_type,sample_once = True)
    return normed_values
    
def loess_graph(density_vectors,called_vectors,data_type = 'DNAse'):
    normed_values = norm_by_graph(density_vectors,called_vectors,
                    method=lowess_norm_pair,sample_number = 25000,data_type=data_type,
                    sample_once = False)
    return normed_values

def tmm_group_norm(density_vectors,called_vectors,data_type = 'DNAse'):
    normed_values =  group_norm(density_vectors,called_vectors,
                        method=tmm_norm_pair,sample_number = 25000,
                        sample_once = True,data_type=data_type)
    return normed_values
    
def loess_group_norm(density_vectors,called_vectors,data_type = 'DNAse'):
    normed_values = group_norm(density_vectors,called_vectors,method=lowess_norm_pair,
                    data_type=data_type,sample_number=25000,sample_once = False)
    return normed_values 

def tmm_geomean_threshold(density_vectors,hard_threshold = 10,filter_on_logratio = (0.25,0.75),
                            pseudocount = 1):
    numbers = np.zeros(len(density_vectors[density_vectors.keys()[0]]))
    gm_means = np.ones(len(numbers))
    for k in density_vectors:
        gm_means *= density_vectors[k]
        numbers += 1
    gm_means = gm_means ** (1/float(len(density_vectors)))
    over_threshold = gm_means > hard_threshold
    
    normed_vectors = {}
    for k in density_vectors:
        peaks2 = density_vectors[k]
        log_ratios = np.log2((gm_means[over_threshold]+pseudocount)/(peaks2[over_threshold]+pseudocount))
        sorted_log_ratios = sorted(log_ratios)

        clow = sorted_log_ratios[int(filter_on_logratio[0]*len(log_ratios))]
        chigh = sorted_log_ratios[int(filter_on_logratio[1]*len(log_ratios))] 
        sampled_pass_filter = np.intersect1d(np.where(log_ratios>clow)[0],
                                    np.where(log_ratios<chigh)[0])

        all_filtered = np.where(over_threshold)[0][sampled_pass_filter]
        weights = (1 / (1/gm_means[all_filtered] + 1/peaks2[all_filtered])) ### should be variance of log(x) assuming poisson based on delta method
        ratio =  2**(np.sum(log_ratios[sampled_pass_filter]*weights)/np.sum(weights))
        normed_vectors[k] = peaks2*ratio
                                                                                
    return normed_vectors 
    
    
def tmm_geomean_norm(density_vectors,called_vectors,sample_number = 30000,data_type = 'DNAse'):
    numbers = np.zeros(len(density_vectors[density_vectors.keys()[0]]))
    log_gm_means = np.zeros(len(numbers))
    for k in density_vectors:
        log_gm_means += np.log(density_vectors[k])
        numbers[called_vectors[k]] += 1
    log_gm_means /= float(len(density_vectors))
    gm_means = np.exp(log_gm_means)
    sorted_indices = np.argsort(gm_means)
    called_twice = np.where(numbers>1)[0]   
    print len(called_twice)
    print np.mean(gm_means)
    sampled = sample_from_range(gm_means[sorted_indices[called_twice]],sorted_indices[called_twice],sample_number = sample_number,
                                percentile_cutoffs = PEAK_THRESHOLDS[data_type],top_cutoff = TOP_THRESHOLDS[data_type])
    normed_vectors = {}
    for k in density_vectors:
        normed_vectors[k] = tmm_norm_pair(gm_means,density_vectors[k],np.ones(len(gm_means)),called_vectors[k],sampled=sampled)
                                                                                
    return normed_vectors    
    
def loess_geomean_norm(density_vectors,called_vectors,sample_number = 30000,data_type = 'DNAse'):
    #### because to the way this is set up, it only works if a
    #### the same general peaks are the same:  look at called twice
    ### might want to change it at some point
    #### also if density vectors are not normalized first, run into problems
    
    numbers = np.zeros(len(density_vectors[density_vectors.keys()[0]]))
    gm_means = np.ones(len(numbers))
    pseudocount = np.min([np.min(density_vectors[d]) for d in density_vectors])
    for k in density_vectors:
        gm_means *= (density_vectors[k]+pseudocount)
        numbers[called_vectors[k]] += 1
    gm_means = gm_means ** (1/float(len(density_vectors)))
    sorted_indexes = np.argsort(gm_means)
    called_twice = np.where(numbers>1)[0]   

    normed_vectors = {}
    for k in density_vectors:
        vector_sort = np.argsort(density_vectors[k][called_twice])
        sampled = sample_from_range(np.sqrt(density_vectors[k][called_twice][vector_sort]),called_twice[vector_sort],sample_number = sample_number,
                                percentile_cutoffs = PEAK_THRESHOLDS[data_type],top_cutoff = TOP_THRESHOLDS[data_type])
        sampled = sampled[gm_means[sampled]!=0]
        uncalled_fraction = np.where(called_vectors[k]==0)[0]
        
        min_gm = np.sqrt(np.min(density_vectors[k][sampled]))
        max_gm = np.sqrt(np.max(density_vectors[k][sampled]))
        number_gm = int(len(sampled)*min_gm/float(max_gm-min_gm))   
        
        if len(uncalled_fraction) <= number_gm:
            uncalled_sample_gm = np.arange(len(uncalled_fraction))
            print('Warning: Few Uncalled Peaks')
        else:
            uncalled_sample_gm = SEED.choice(len(uncalled_fraction),number_gm,replace = False) 
            update_seed()
        sample_subset = np.append(sampled,uncalled_fraction[uncalled_sample_gm])

        pseudocount = min(pseudocount,np.min(density_vectors[k][np.where(density_vectors[k]>0)[0]]))      
        new1 = sm.nonparametric.lowess(gm_means[sample_subset],density_vectors[k][sample_subset],
                                    return_sorted=False,it=4,frac = 0.4,delta = 0.01*max(density_vectors[k][sample_subset]))        
                             
        interpolated = return_interpolated(density_vectors[k],new1,sample_subset,pseudocount)
        #import matplotlib.pyplot as plt
        #plt.figure()
        #plt.plot(density_vectors[k][sample_subset],gm_means[sample_subset],marker='o',linestyle='',alpha = 0.2)
        #plt.plot(sorted(density_vectors[k][sample_subset]),sorted(new1),color='r')
        #import pdb
        #pdb.set_trace()
        normed_vectors[k] = interpolated
                                                                                
    return normed_vectors
    
def group_norm(density_vectors,called_vectors,method=None,data_type = 'DNAse',
                sample_number = 10000,sample_once = False):  
    sorted_indexes,gm_means,numbers = return_ordered_indexes(density_vectors,called_vectors)
    called_twice = np.where(numbers>1)[0]    
    
    if sample_once:
        sampled = sample_from_range(gm_means[sorted_indexes[called_twice]],sorted_indexes[called_twice],sample_number = sample_number,
                                percentile_cutoffs = PEAK_THRESHOLDS[data_type],top_cutoff = TOP_THRESHOLDS[data_type])
    else:    
        sampled = None
    
    pseudocount = 1.
    for k in density_vectors:
        pseudocount = min(pseudocount,np.min(density_vectors[k][np.where(density_vectors[k]>0)[0]]))      
    
    ref_distances = {}
    for r in density_vectors:
        ref_distances[r] = np.mean((density_vectors[r][called_vectors[r]]-gm_means[called_vectors[r]])**2)
    ref_exp = sorted([r for r in ref_distances],key = lambda k:ref_distances[r])[0]
    
    normed_vectors = {}
    for k in density_vectors:
        if k == ref_exp:
            normed_vectors[k] = density_vectors[k]
        else:
            if FILTER_CALLED[data_type] and sampled != None:
                current_sample = filter_on_peaks(sampled,called_vectors[ref_exp],called_vectors[k])
            else:
                current_sample = sampled          
            normed_vectors[k] = method(density_vectors[ref_exp],density_vectors[k],called_vectors[ref_exp],called_vectors[k],sampled=current_sample,
                                        filter_on_logratio = TMM_RATIOS[data_type],pseudocount = pseudocount)   
    return normed_vectors

def norm_by_graph(density_vectors,called_vectors,method=None,
                        sample_number = 20000,data_type = 'DNAse',sample_once = False):
    sorted_indexes,gm_means,numbers = return_ordered_indexes(density_vectors,called_vectors)
    called_twice = np.where(numbers>1)[0]    
    
    if sample_once:
        sampled = sample_from_range(gm_means[sorted_indexes[called_twice]],sorted_indexes[called_twice],sample_number = sample_number,
                                percentile_cutoffs = PEAK_THRESHOLDS[data_type],top_cutoff = TOP_THRESHOLDS[data_type])
    else:    
        sampled = None
    
    pseudocount = 1.
    for k in density_vectors:
        pseudocount = min(pseudocount,np.min(density_vectors[k][np.where(density_vectors[k]>0)[0]]))      
    
    ref_distances = {}
    for r in density_vectors:
        ref_distances[r] = np.mean((density_vectors[r][called_vectors[r]]-gm_means[called_vectors[r]])**2)
    ref_exp = sorted([r for r in ref_distances],key = lambda k:ref_distances[r])[0]
                            
    all_means = sorted([(np.mean(density_vectors[k]),k) for k in density_vectors])
    reference_exp = all_means[-1][1]    
    normalization_values  = {reference_exp:density_vectors[reference_exp]}       
    for exp1 in density_vectors:
        if exp1 == reference_exp:
            continue
        path_values = np.zeros(len(density_vectors[reference_exp]))
        total_weight = np.zeros(len(path_values))
        for exp2 in density_vectors:
            if exp2 == reference_exp:
                continue
            key1 = (exp1,exp2)
            if exp2 == exp1:
                if FILTER_CALLED[data_type]:
                    current_sampled = filter_on_peaks(sampled,called_vectors[reference_exp],called_vectors[exp1])
                else:
                    current_sampled = sampled
                new_vector = method(density_vectors[reference_exp],density_vectors[exp1],called_vectors[reference_exp],called_vectors[exp1],
                                    current_sampled,filter_on_logratio = TMM_RATIOS[data_type],
                                    pseudocount = pseudocount)
                                        
                lratio = np.log(new_vector/density_vectors[exp1])                
                weight = 1/np.exp(np.mean(lratio[np.isfinite(lratio)]))**2
                path_values += new_vector*weight
                total_weight += weight
            
            else:
                if FILTER_CALLED[data_type]:
                    current_sampled = filter_on_peaks(sampled,called_vectors[exp1],called_vectors[exp2])  
                else:
                    current_sampled = sampled
                intermed = method(density_vectors[exp2],density_vectors[exp1],called_vectors[exp2],called_vectors[exp1],
                                    current_sampled,filter_on_logratio = TMM_RATIOS[data_type],
                                    pseudocount = pseudocount)                
                lratio1 = np.log(intermed/density_vectors[exp1])
                if FILTER_CALLED[data_type]:
                    current_sampled = filter_on_peaks(sampled,called_vectors[exp2],called_vectors[reference_exp])     
                else:
                    current_sampled = sampled
                final = method(density_vectors[reference_exp],intermed,called_vectors[reference_exp],called_vectors[exp1],
                                current_sampled,
                                filter_on_logratio = TMM_RATIOS[data_type],
                                pseudocount = pseudocount)                                
                lratio2 = np.log(final/intermed)
                weight = 1/(np.exp(np.mean(lratio1[np.isfinite(lratio1)]))**2 + 
                            np.exp(np.mean(lratio2[np.isfinite(lratio2)]))**2)
                path_values += final*weight
                total_weight += weight
        normalization_values[exp1] = path_values/total_weight
    return normalization_values
    
    
def sample_from_range(sorted_geometric_means,sorted_indexes,percentile_cutoffs = [.05,1],
                        sample_number=10000,tolerance_percentage=0.01,subset = None,top_cutoff = 0.05,hard_cutoffs = None):
    cut_indexes = [int((len(sorted_geometric_means)-1)*r) for r in percentile_cutoffs]
    while cut_indexes[-1] >= (len(sorted_geometric_means)-1):
        cut_indexes[-1] -= 1
    cut_indexes[-1] = filter_by_range(sorted_geometric_means[0:cut_indexes[-1]],top_cutoff)  

    if hard_cutoffs != None:
        cut_indexes[0] = np.searchsorted(sorted_geometric_means,hard_cutoffs[0])
        cut_indexes[1] = np.searchsorted(sorted_geometric_means,hard_cutoffs[1])
    sorted_log_heights = sorted_geometric_means
    cutoffs = [sorted_log_heights[i] for i in cut_indexes]

    tolerance = (cutoffs[1] - cutoffs[0])* tolerance_percentage

    sample_slices = []
    indexes = set([])
    min_mean = sorted_log_heights[0]
    max_mean = sorted_log_heights[-1]    
    for r in xrange(sample_number):
        hcenter = SEED.random_sample()*(cutoffs[1]-cutoffs[0]) + cutoffs[0]
        update_seed()
        #random.uniform(cutoffs[0],cutoffs[1])
        sample_slices.append((max(min_mean,hcenter-tolerance),min(max_mean,hcenter+tolerance)))
    tolerable_min_index = np.searchsorted(sorted_log_heights,[t[0] for t in sample_slices],side = 'left')
    tolerable_max_index = np.searchsorted(sorted_log_heights,[t[1] for t in sample_slices],side = 'right')-1

    
    for i,x in enumerate(tolerable_min_index):
        if x > tolerable_max_index[i]:
            continue
        if subset != None:
            acceptable = np.intersect1d(np.arange(x,tolerable_max_index[i]+1),subset)
            if len(acceptable) == 0:
                continue
            indexes.add(int(SEED.choice(acceptable)))
            update_seed()
        else:
            acceptable = np.setdiff1d(np.arange(x,tolerable_max_index[i]+1),indexes)
            if len(acceptable) == 0:
                continue
                #print x,tolerable_max_index[i]
            indexes.add(int(SEED.choice(acceptable)))
            update_seed()
    index_random = sorted_indexes[list(indexes)]
    return index_random

def return_ordered_indexes(density_vectors,called_vectors):
    counts = np.zeros(len(density_vectors[density_vectors.keys()[0]]))
    gm_means = np.ones(len(counts))
    for k in density_vectors:
        gm_means[called_vectors[k]] *= density_vectors[k][called_vectors[k]]
        counts[called_vectors[k]] += 1
    gm_means = gm_means ** (1/counts)
    sorted_indexes = np.argsort(gm_means)
    return sorted_indexes,gm_means,counts
######
def filter_on_peaks(core_subset,called1,called2):###    
    good = np.intersect1d(np.where(called1>0)[0],np.where(called2>0)[0])
    return np.intersect1d(core_subset,good)
    
def filter_by_range(sorted_list,threshold):
    LIMIT = 40000
    from scipy.stats import expon
    if int(len(sorted_list)*.99) > LIMIT:
        subset = SEED.choice(int(.99*len(sorted_list)),LIMIT)
        update_seed()
    else:
        subset = np.arange(int(len(sorted_list)*.99))
    fitted = expon.fit(sorted_list[0:int(.99*len(sorted_list))][subset])
    r = -1
    while len(sorted_list)+r > 0:
        pv = expon.pdf(sorted_list[r],fitted[0],fitted[1])*len(sorted_list)
        if pv > threshold:
            break
        r -= 1
    return len(sorted_list) + r    
     
    
def tmm_norm_pair(peaks1,peaks2,called1,called2,sampled = None,pseudocount = 1.,filter_on_logratio = (0.25,0.75)):

    if sampled == None:
        sort_order = np.argsort(peaks1*peaks2)
        sampled = sample_from_range(np.sqrt(peaks1*peaks2)[sort_order],sort_order)
    ##### need to add in log ratio cutoffs
    log_ratios = np.log2((peaks1[sampled]+pseudocount)/(peaks2[sampled]+pseudocount))
    sorted_log_ratios = sorted(log_ratios)
    if filter_on_logratio:
        clow = sorted_log_ratios[int(filter_on_logratio[0]*len(log_ratios))]
        chigh = sorted_log_ratios[int(filter_on_logratio[1]*len(log_ratios))] 
        sampled_pass_filter = np.intersect1d(np.where(log_ratios>clow)[0],
                                    np.where(log_ratios<chigh)[0])
    else:
        sampled_pass_filter = np.where(np.isfinite(log_ratios))[0]
    all_filtered = sampled[sampled_pass_filter]
    weights = (1 / (1/peaks1[all_filtered] + 1/peaks2[all_filtered])) ### should be variance of log(x) assuming poisson based on delta method
    ratio =  2**(np.sum(log_ratios[sampled_pass_filter]*weights)/np.sum(weights))
    return peaks2*ratio
   

def lowess_norm_pair(peaks1,peaks2,called1,called2,sampled = None,pseudocount = 1.,
                        filter_on_logratio = False):
    geomean = np.sqrt(peaks1*peaks2)
    if sampled == None:
        called_both = np.where(called1*called2 > 0)[0]
        peaks1_s = np.argsort(peaks1[called_both])
        peaks2_s = np.argsort(peaks2[called_both])
        sampled_1 = sample_from_range(np.sqrt(peaks1[called_both][peaks1_s]),
                                    peaks1_s,percentile_cutoffs=[0,1],sample_number=15000,top_cutoff=0.025)
        sampled_2 = sample_from_range(np.sqrt(peaks2[called_both][peaks2_s]),
                                    peaks2_s,percentile_cutoffs=[0,1],sample_number=15000,top_cutoff=0.025)               
    
    uncalled_fraction = np.where(called1*called2 == 0)[0]
        
    sorted_uncalled_1 = np.argsort(peaks1[uncalled_fraction])
    min_1 = np.min(peaks1[called_both][sampled_1])
    max_1 = np.max(peaks1[called_both][sampled_1])
    number_1 = int(len(sampled_1)*min_1/float(max_1-min_1))
    sorted_uncalled_2 = np.argsort(peaks2[uncalled_fraction])
    min_2 = np.min(peaks2[called_both][sampled_2])
    max_2 = np.max(peaks2[called_both][sampled_2])
    number_2 = int(len(sampled_2)*min_2/float(max_2-min_2))    
    
    #uncalled_sample1 = sample_from_range(peaks1[uncalled_fraction][sorted_uncalled_1],sorted_uncalled_1,
    #                                sample_number=number_1,hard_cutoffs=(0,min_1)) 
    #uncalled_sample2 = sample_from_range(peaks2[uncalled_fraction][sorted_uncalled_2],sorted_uncalled_2,
    #                                sample_number=number_2,hard_cutoffs=(0,min_2))  
    if len(uncalled_fraction) <= number_1:
        uncalled_sample1 = np.arange(len(uncalled_fraction))
        print 'Warning: Few Uncalled Peaks'
    else:
        uncalled_sample1 = SEED.choice(len(uncalled_fraction),number_1,replace = False) 
        update_seed()
    if len(uncalled_fraction) <= number_2:
        print 'Warning: Few Uncalled Peaks'
        uncalled_sample2 = np.arange(len(uncalled_fraction))
    else:
        uncalled_sample2 = SEED.choice(len(uncalled_fraction),number_2,replace = False)    
        update_seed()
    #### need to add in all the way down to zero
    sampled_1 = np.append(called_both[sampled_1],uncalled_fraction[uncalled_sample1])
    sampled_2 = np.append(called_both[sampled_2],uncalled_fraction[uncalled_sample2])
    
    #import matplotlib
    #matplotlib.use('Agg')
    #import matplotlib.pyplot as plt    
  
    #new0 = sm.nonparametric.lowess(peaks1[all_todo],peaks2[all_todo],
    #                                return_sorted=False,it=4,frac = 0.5,delta = 0.01*max(peaks2[all_todo]))
    #new1 = sm.nonparametric.lowess(geomean[all_todo],peaks1[all_todo],
    #                                return_sorted=False,it=4,frac = 0.5,delta = 0.01*max(peaks1[all_todo]))        
    #new2 = sm.nonparametric.lowess(geomean[all_todo],peaks2[all_todo],
    #                                return_sorted=False,it=4,frac = 0.5,delta = 0.01*max(peaks2[all_todo]))  
    new1 = sm.nonparametric.lowess(geomean[sampled_1],peaks1[sampled_1],
                                    return_sorted=False,it=4,frac = 0.4,delta = 0.01*max(peaks1[sampled_1]))        
    new2 = sm.nonparametric.lowess(geomean[sampled_2],peaks2[sampled_2],
                                    return_sorted=False,it=4,frac = 0.4,delta = 0.01*max(peaks2[sampled_2]))
                                    
   
    #interpolated1 = return_interpolated(peaks1,new1,sampled_1,pseudocount)
    interpolated2 = return_interpolated(peaks2,new2,sampled_2,pseudocount)
    
    sampled_order = np.argsort(new1)
    interpolated3 =  np.interp(interpolated2,new1[sampled_order],peaks1[sampled_1][sampled_order])    
    interpolated3 = extrapolate(interpolated3,interpolated2,new1,peaks1[sampled_1])
    
    return interpolated3

    
def quantile_normalization(density_vectors):
    sorted_args = {}
    for d in density_vectors:
        sorted_args[d] = np.argsort(density_vectors[d])
    mean_vector = np.zeros(len(density_vectors[d]))
    for d in sorted_args:
        mean_vector += density_vectors[d][sorted_args[d]]
    mean_vector /= len(density_vectors.keys())
    normed_vectors = {}
    for d in density_vectors:
        normed_vectors[d] = np.zeros(len(density_vectors[d]))
        normed_vectors[d][sorted_args[d]] = mean_vector
    return normed_vectors    
    
def extrapolate(interpolated,to_extrap,base,predict):
    EXTRAPOLATION_FRACTION = 0.1
    sample_order = np.argsort(base)
    min_index_range = (sample_order[0],sample_order[int(EXTRAPOLATION_FRACTION*len(sample_order))])
    # min_index_range = (sample_order[base[sample_order != 0]][0], sample_order[int(EXTRAPOLATION_FRACTION * len(sample_order))])
    max_index_range = (sample_order[-1],sample_order[-1*int(EXTRAPOLATION_FRACTION*len(sample_order))])
    slope_min = (predict[min_index_range[0]])/float(base[min_index_range[0]])
    slope_max = (predict[max_index_range[0]]-predict[max_index_range[1]])/float(base[max_index_range[0]]-base[max_index_range[1]])

    under_min = np.where(to_extrap < base[min_index_range[0]])[0]
    over_max = np.where(to_extrap > base[max_index_range[0]])[0]
    under_lambda = lambda r: (r - base[min_index_range[0]])*slope_min+predict[min_index_range[0]]
    over_lambda = lambda r: (r - base[max_index_range[0]])*slope_max+predict[max_index_range[0]]
    interpolated[under_min] = under_lambda(to_extrap[under_min])
    interpolated[over_max] = over_lambda(to_extrap[over_max])
    interpolated[np.where(interpolated < 0.0)[0]] = 0.0

    return interpolated
    
def return_interpolated(peaks,lowess_est,sample,pseudocount = 1):
    unique,unique_indices = np.unique(peaks[sample],return_index = True)
    interpolated =  np.interp(peaks,peaks[sample][unique_indices],lowess_est[unique_indices])  
    interpolated = extrapolate(interpolated,peaks,peaks[sample],lowess_est)
    return interpolated
    min_index = sampled_order[0]
    max_index = sampled_order[-1]
    min_value_ratio = (lowess_est[min_index]+pseudocount)/(peaks[sample][min_index]+pseudocount) 
    max_value_ratio = lowess_est[max_index]/peaks[sample][max_index]  ### if zero giving issues
    over_max = np.where(peaks>peaks[sample][max_index])[0]
    interpolated[over_max] = max_value_ratio * peaks[over_max]
    under_min = np.where(peaks<peaks[sample][min_index])[0]
    interpolated[under_min] = np.maximum(0,min_value_ratio * peaks[under_min])
    interpolated[np.where(interpolated<0)[0]] = 0.0  
    return interpolated
   

######################################## code junkyard past here    
def group_norm_geometric(density_vectors,called_vectors,method=None,data_type = 'DNAse',
                sample_number = 10000):   
    #### this doesn't really work well. Hard to do will varying #s of times called
    sorted_indexes,gm_means,numbers = return_ordered_indexes(density_vectors,called_vectors)
    called_twice = np.where(numbers>1)[0]
    sampled = sample_from_range(gm_means[sorted_indexes[called_twice]],sorted_indexes[called_twice],sample_number = sample_number,
                                percentile_cutoffs = PEAK_THRESHOLDS[data_type],top_cutoff = TOP_THRESHOLDS[data_type])
    print len(sampled)
    print len(called_twice)
    pseudocount = 1.
    for k in density_vectors:
        pseudocount = min(pseudocount,np.min(density_vectors[k][np.where(density_vectors[k]>0)[0]]))      
    
    normed_vectors = {}
    for k in density_vectors:
        if FILTER_CALLED[data_type]:
            current_sample = np.intersect1d(sampled,np.where(called_vectors[k]>0)[0])
        else:
            current_sample = sampled          
        print len(current_sample)
        normed_vectors[k] = method(gm_means,density_vectors[k],current_sample,
                                        filter_on_logratio = TMM_RATIOS[data_type],pseudocount = pseudocount)   
    return normed_vectors 


class Fake:
    def __init__(self):
        pass
    
def run_norm(names,peak_list,tag_list,location = '/home/jlazar/tmp/master_list.bed'):
    import refList
    import gridmap
    fakelist = []
    for i,n in enumerate(names):
        f = Fake()
        f.name = n
        f.peaks_filepath = peak_list[i]
        f.tag_filepath = tag_list[i]
        fakelist.append(f)
    rL = refList.Master_List(fakelist, location)
    height_jobs = []
    shared_jobs = []
    ziptemps = []
    for i,n in enumerate(names):
        ziptemp = bedtools.make_tmpfile()
        ziptemps.append(ziptemp)
        bedtools.uncompress(peak_list[i], ziptemp)
        height_jobs.append(gridmap.Job(bedtools.bedmap, ['tempfile', '--max', location, tag_list[i]]))
        shared_jobs.append(gridmap.Job(bedtools.bedmap, ['tempfile', '--indicator', location, ziptemp]))
    height_out = gridmap.process_jobs(height_jobs)
    shared_out = gridmap.process_jobs(shared_jobs)
    weighted_vectors = {}
    shared_vectors = {}
    for i,n in enumerate(names):
        weighted_vectors[n] = np.loadtxt(height_out[i])
        shared_vectors[n] = np.loadtxt(shared_out[i])
    for z in ziptemps:
        os.remove(z)
    return weighted_vectors,shared_vectors
    return tpp_norm_group(weighted_vectors,shared_vectors)
    
def trim_vector(vector,low_p,high_p):
    sorted_vec = sorted(vector)
    low_ind = int(low_p*len(sorted_vec))
    high_ind = int(high_p*len(sorted_vec))
    low = sorted_vec[low_ind]
    high = sorted_vec[high_ind]
    good_ind = np.intersect1d(np.where(vector>low)[0],np.where(vector<high)[0])
    return vector[good_ind]

def split2percentiles(vector,number=100):
    to_sort = range(len(vector))
    to_sort.sort(key=lambda k:vector[k])
    each_section = len(to_sort)/float(number)
    quartile_indices = []
    start = 0
    for r in range(number):
        next_ind = int(round(start+each_section))
        quartile_indices.append([c for c in to_sort[int(round(start)):next_ind]])
        start = start+each_section
    return quartile_indices    

def interpolate(index_dict,vector):
    ### create a lambda
    sorted_x = np.array(sorted(index_dict.keys()))
    sorted_y = np.array([index_dict[k] for k in sorted_x])
    tmp =  np.interp(vector,sorted_x,sorted_y,right = np.nan)
    extrapolate = np.mean(sorted_y[-10:]/sorted_x[-10:])
    tmp[np.isnan(tmp)] *= extrapolate
    #this really needs to be done
    
def index_vector(vector):
    all_values = set(vector)
    index ={}
    for v in all_values:
        index[v] = np.where(vector==v)[0]
    return index
def get_outliers(index_dict,total_length):
    passed = 0
    outlier_threshold = (1-0.0001)*total_length
    sorted_keys = sorted(index_dict.keys())
    for i,k in enumerate(sorted_keys):
        passed += len(index_dict[k])
        if passed > outlier_threshold:
            break
    return set(sorted_keys[i+1:])    

def tpp_norm_old(peaks1,peaks2,called1,called2,sort_order):
    quantiles1 = split2percentiles(peaks1)
    shared_percents = [np.mean(shared[q]) for q in quantiles1]
    min_percentile = .25
    max_sp = max(shared_percents)
    for i,sp in enumerate(shared_percents):
        if max_sp-sp<0.05:
            min_percentile = i/float(len(shared_percents))
            break
    max_percentile = max(.95,min_percentile+0.01)
    
    sorted_1 = sorted(peaks1)
    sorted_2 = sorted(peaks2)
    cutoffs_1 = (sorted_1[int(min_percentile*len(sorted_1))],sorted_1[int(max_percentile*len(sorted_1))])
    cutoffs_2 = (sorted_2[int(min_percentile*len(sorted_2))],sorted_2[int(max_percentile*len(sorted_2))])
    p1 = np.intersect1d(np.where(peaks1>cutoffs_1[0])[0],np.where(peaks1<cutoffs_1[1])[0])
    p2 = np.intersect1d(np.where(peaks2>cutoffs_2[0])[0],np.where(peaks2<cutoffs_2[1])[0])
    all_filtered = np.intersect1d(p1,p2)
    
    log_ratios = np.log2(peaks1[all_filtered]/peaks2[all_filtered])
    weights = 1/peaks1[all_filtered] + 1/peaks2[all_filtered]
    return 2**(sum(log_ratios*weights)/sum(weights))    
    
def find_intercept(mean_vector,data_vector,subset):
    cutoff_adjust = 0.00001 # this is a small factor to fix bug if too many of same values
    subset = np.where(subset>0)[0]
    log2s = np.log2(data_vector[subset]/mean_vector[subset])
    sorted_ls = sorted(log2s)
    cutoffs = (sorted_ls[int(.25*len(log2s))]-cutoff_adjust,
               sorted_ls[int(.75*len(log2s))]+cutoff_adjust)
    ns = np.intersect1d(np.where(log2s>cutoffs[0])[0],np.where(log2s<cutoffs[1])[0])

    model = sm.RLM(mean_vector[subset][ns],sm.add_constant(data_vector[subset][ns]),M=sm.robust.norms.HuberT())
    results = model.fit()

    return results.params[0]

def set_by_pairs(all_exps,pairwise_dict):
    normalization_factors = {}
    reference_exp = all_exps[0]
    normalization_factors[reference_exp] = 1
    for exp1 in all_exps:
        if exp1==reference_exp:
            continue
        path_values = 0
        for exp2 in all_exps:
            if exp2 == reference_exp:
                continue
            key1 = tuple(sorted([reference_exp,exp2],key=lambda k:all_exps.index(k)))
            if exp2 == exp1:
                path_values += (2*pairwise_dict[key1][exp1][1]/pairwise_dict[key1][reference_exp][1])
            else:
                key2 = tuple(sorted([exp1,exp2],key=lambda k:all_exps.index(k)))
                path_values += ((pairwise_dict[key1][exp2][1]/pairwise_dict[key1][reference_exp][1])*
                                (pairwise_dict[key2][exp1][1]/pairwise_dict[key2][exp2][1]))
        normalization_factors[exp1] = path_values/float(len(all_exps))
    return normalization_factors
       
    
def tpp_norm_group(weighted_vector_dict,shared_vector_dict):
    exp_number = len(weighted_vector_dict.keys())
    heights_vector = np.zeros(len(weighted_vector_dict[weighted_vector_dict.keys()[0]]))+1
    count_vector = np.zeros(len(heights_vector))
    for exp in weighted_vector_dict:
        called_peaks = np.where(shared_vector_dict[exp]>0)[0]
        #heights_vector[called_peaks] *= weighted_vector_dict[exp][called_peaks]
        heights_vector *= weighted_vector_dict[exp]
        count_vector[called_peaks] += 1.0
    mean_vector = heights_vector**(1/count_vector)
    real_peaks = np.where(count_vector>0)[0]
    mean_vector = mean_vector[real_peaks]
    count_vector = count_vector[real_peaks]

    quantiles = split2percentiles(mean_vector)
    shared_percents = [np.mean(count_vector[q]/exp_number) for q in quantiles]
    means_percents = [np.mean(mean_vector[q]) for q in quantiles]
    min_percentile = .25
    max_sp = max(shared_percents)

    for i,sp in enumerate(shared_percents):
        if max_sp-sp<0.1:
            min_percentile = i/float(len(shared_percents))
            break
   # print max_sp,min_percentile,shared_percents
    max_percentile = max(shared_percents.index(max_sp)/100.,min_percentile+0.01)

    sorted_ref = sorted(mean_vector)
    cutoffs = (sorted_ref[int(min_percentile*len(sorted_ref))],sorted_ref[int(max_percentile*len(sorted_ref))])
    height_filtered = np.intersect1d(np.where(mean_vector>cutoffs[0])[0],np.where(mean_vector<cutoffs[1])[0])
    norm_factors = {}
    for exp in weighted_vector_dict:
        intercept = 0#find_intercept(mean_vector,weighted_vector_dict[exp],shared_vector_dict[exp])
        all_filtered = np.intersect1d(height_filtered,np.where(shared_vector_dict[exp]>0)[0])
        log_ratios = np.log2((weighted_vector_dict[exp][all_filtered]-intercept)/mean_vector[all_filtered])
        weights = (np.sum(mean_vector)-mean_vector[all_filtered])/(mean_vector[all_filtered]*np.sum(mean_vector))
        sorted_ratios = sorted(log_ratios)
        ctop = sorted_ratios[int(TMM_RATIO_TOP*len(sorted_ratios))]
        cbottom = sorted_ratios[int(TMM_RATIO_BOTTOM*len(sorted_ratios))]
        final_trim = np.intersect1d(np.where(log_ratios>cbottom)[0],np.where(log_ratios<ctop)[0])
        norm_factors[exp] =  (intercept,1/(2**(sum(log_ratios[final_trim]*weights[final_trim])/sum(weights[final_trim]))))      
    return norm_factors    
    
def create_normalization_file(normed_vectors,raw_vectors,name_order,normalization_file = '',
                                header = False):      
    pseudocount_raw = 1. #### this is not robust,, right now dealing with small cp/n sucks
    size_matrix = []
    for i in name_order:
        n = normed_vectors[i]
        cp = raw_vectors[i]  
        #median_index = np.median
        #median_ratio = np.median(cp>)/np.median(n)
        median_ratio = np.median(cp[np.logical_and(n>0,cp>0)]/n[np.logical_and(n>0,cp>0)])
        norm_pseudocount = float(np.min(n[np.logical_and(n>0,cp>0)]))
        raw_pseudocount = float(np.min(cp[np.logical_and(n>0,cp>0)]))
        scales = cp/n
        potential_issue = np.logical_or(n<1,cp<1)
        scales[potential_issue] = median_ratio ## ? median ratio
        size_matrix.append(scales)

    size_matrix = np.array(size_matrix)
    size_matrix[np.isnan(size_matrix)]=1.0
    size_matrix[np.logical_not(np.isfinite(size_matrix))]=1.0

    ### geometric mean
    #geom_means = np.exp(np.mean(np.log(size_matrix),axis=0))
    geom_means = np.exp(np.mean(np.log(size_matrix)))
    size_matrix /= geom_means
    if normalization_file == '':
        out_file = bedtools.make_tmpfile()
    else:
        out_file = normalization_file
    if not header:
        np.savetxt(out_file,np.transpose(size_matrix),'%.5e',delimiter = '\t')
    else:
        np.savetxt(out_file,np.transpose(size_matrix),'%.5e',header = '\t'.join(name_order),comments = '',
                        delimiter = '\t')
    return out_file    