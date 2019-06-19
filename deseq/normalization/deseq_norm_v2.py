import argparse
import logging
import numpy as np
import os
import time,tempfile
from scipy.stats import expon,spearmanr

import smoothers_lowess as sm

# script log created as deseq2_{today}.log in working directory
today = time.strftime("%d%m%Y")
FORMAT = "[%(name)s | Line:%(lineno)s| Function:%(funcName)s ] %(message)s"
logging.basicConfig(filename="deseq2_{}.log".format(today), level=logging.INFO, filemode='w',
                    format=FORMAT)
log = logging.getLogger('normalization_script')


SEED_NUMBER = 1832245
SEED = np.random.RandomState(SEED_NUMBER)


PEAK_OUTLIER_THRESHOLD = 0.999
DELTA_FRACTION = 0.001
CORRELATION_LIMIT = 0.8
CV_FRACTION = 0.33
SAMPLE_NUMBER = 75000 
BIN_NUMBER = 100


def return_averages(density_vectors,called_vectors):
    counts = np.zeros(len(density_vectors[density_vectors.keys()[0]]))
    gm_means = np.ones(len(counts))
    for k in density_vectors:
        gm_means[called_vectors[k]] += np.log(density_vectors[k][called_vectors[k]])
        counts[called_vectors[k]] += 1
    #gm_means = np.exp(gm_means * float(1/counts))
    
    return gm_means,counts

def outlier_limit(array_values,threshold,limit = SAMPLE_NUMBER):
    if len(array_values) > limit:
        subset = SEED.choice(len(array_values),limit)
    else:
        subset = np.arange(len(array_values))
    fitted = expon.fit(array_values[subset])
    return expon(fitted[0],fitted[1]).ppf(1-(1-threshold)/len(array_values))
    
    r = -1
    while len(sorted_list)+r > 0:
        pv = expon.pdf(sorted_list[r],fitted[0],fitted[1])*len(sorted_list)
        if pv > threshold:
            break
        r -= 1
       
    return len(sorted_list) + r

def select_uniform(vector_value,number = SAMPLE_NUMBER,bins = BIN_NUMBER,top_percentage = PEAK_OUTLIER_THRESHOLD,log = False,random = False,
                    verboten = None,sample_method = 'raw'):
    if verboten is None: verboten = np.array([])
    not_verboten = np.setdiff1d(np.arange(len(vector_value)),verboten)
    
    if sample_method == 'log':
        vls = np.log(vector_value[not_verboten]+1)
        good = not_verboten
    else:
        max_value = outlier_limit(vector_value,top_percentage)
        good = np.where(vector_value<max_value)[0]
        good = np.intersect1d(good,not_verboten)
        vls = vector_value[good]
    
    if sample_method == 'random':
        return SEED.choice(len(vector_value),size = min(len(vector_value),number),replace = False).astype(int)    
    
    bin_size = (np.max(vls)-np.min(vls))/float(bins)
    number_each = int(np.ceil(number/bins))
    final_index = np.array([])
    
    for i in np.arange(np.min(vls),np.max(vls),bin_size):
        window_min = i
        window_max = i+bin_size
        
        good_ind = np.where((vls>=window_min) & (vls<window_max))[0]
        if len(good_ind) == 0: continue
        final_index = np.union1d(final_index,SEED.choice(good_ind,size = min(len(good_ind),number_each),replace = False))

    return good[final_index.astype(int)]
    
    
def lowess_group_norm(density_vectors, called_vectors,
               sample_number=SAMPLE_NUMBER, sample_once=False,sample_method = 'raw'):
    
    gm_means, numbers = return_averages(density_vectors, called_vectors)
    decent_peaks = get_peak_subset(gm_means,numbers)

    if sample_once:
        sampled = select_uniform(gm_means[decent_peaks],number = sample_number,top_percentage = PEAK_THRESHOLD,
                            sample_method = sample_method)
        sampled = decent_peaks[sampled]
    else:
        sampled = None

    pseudocounts = {d: np.min(density_vectors[d][density_vectors[d]>0]) for d in density_vectors}

    ref_distances = {}
    for r in density_vectors:
        ref_distances[r] = np.mean((density_vectors[r] - gm_means) ** 2)
    ref_exp = sorted([r for r in ref_distances], key=lambda k: ref_distances[r])[0]

    normed_vectors = {}
    for k in density_vectors:
        if k == ref_exp:
            normed_vectors[k] = density_vectors[k]
        else:
            if not_sample_once:
                decent_peaks = np.where((called_vectors[k] + called_vectors[ref_exp])>1)[0]
                gm_means = np.exp(np.log(density_vectors[k])/2.+np.log(density_vectors[ref_exp])/2.)
                sampled = select_uniform(gm_means[decent_peaks],number = sample_number,top_percentage = PEAK_THRESHOLD,
                            sample_method = sample_method)
                            
            normed_vectors[k] = lowess_norm_pair(density_vectors[ref_exp], density_vectors[k], called_vectors[ref_exp],
                                       called_vectors[k], sampled = sampled,
                                       pseudocount=pseudocounts)
    return normed_vectors

### sets the values out of the range to the edge values.  For 
### logged data, equivalent to linear approximation between 0 and the point
def extrapolate(interpolated,to_extrap,base,predict):

    sample_order = np.argsort(base)
    min_value = np.min(base)
    max_value = np.max(base)
    
    under_min = np.where(to_extrap < min_value)[0]
    over_max = np.where(to_extrap > max_value)[0]
    
    min_predict = predict[sample_order[0]]
    max_predict = predict[sample_order[-1]]
    
    interpolated[under_min] = min_predict
    interpolated[over_max] = max_predict

    return interpolated

def return_interpolated(peaks,lowess_est,sample):
    unique,unique_indices = np.unique(peaks[sample],return_index = True)
    interpolated =  np.interp(peaks,peaks[sample][unique_indices],lowess_est[unique_indices])
    interpolated = extrapolate(interpolated,peaks,peaks[sample],lowess_est)
    return interpolated

def lowess_norm_pair(peaks1, peaks2, called1, called2, sampled=None, pseudocount=1,sample_method = 'random'):
    log_geomean = (np.log(peaks1)+np.log(peaks2))/2.
    if type(pseudocount) is list:
        dif = np.log(peaks2+pseudocount[1])-np.log(peaks1+pseudocount[1])
        pseudocount = pseudocount[1]
    else:
        dif = np.log(peaks2+pseudocount)-np.log(peaks1+pseudocount)
    
    if sampled == None:
        sampled = select_uniform(np.exp(log_geomean),number = SAMPLE_NUMBER,bins = BIN_NUMBER,top_percentage = PEAK_OUTLIER_THRESHOLD,sample_method = sample_method)

    
    sampled = sampled[(peaks1[sampled] != 0) & (peaks2[sampled] != 0)]
    xvalues = np.log(peaks1)
    cv_fraction = choose_fraction_cv(xvalues,dif,sampled)
    new_values = sm.lowess(dif[sampled],xvalues[sampled],
                                   return_sorted=False, it=4, frac = cv_fraction)
    interpolated2 = return_interpolated(xvalues, new_values, sampled)
    return peaks2/np.exp(interpolated2)

### choose smoothing parameter for loess by cross-validation
def choose_fraction_cv(xdata,ydata,sample,start = 0.1,end = 0.8,step = 0.1,delta = None):
    in_range = np.where((xdata>=np.min(xdata[sample])) & (xdata<=np.max(xdata[sample])))[0]  ## only xvalidate in range
    cv_sample = np.setdiff1d(in_range,sample)
    if len(cv_sample) < CV_FRACTION*len(sample)/float(1-CV_FRACTION):
        cv_number = int(len(in_range) * CV_FRACTION)
        cv_sample = SEED.choice(in_range,cv_number,replace = False)
        sample = np.setdiff1d(in_range,cv_sample)
        
    
    min_err = np.inf
    best_frac = 0
    if delta == None:
        delta = DELTA_FRACTION * np.percentile(xdata,99)
    
    for frac in np.arange(start,end+step,step):
        new_values = sm.lowess(ydata[sample],xdata[sample],
                                   return_sorted=False, it=4, frac = frac,delta = delta)
        interpolated2 = return_interpolated(xdata, new_values, sample)
        if np.isnan(np.max(interpolated2)):
            err = np.inf
        else:
            err = np.mean((interpolated2[cv_sample]-ydata[cv_sample])**2)
        if err < min_err:
            min_err = err
            best_frac = frac
    return best_frac
        

### normalize all peaks to the geometric mean
### this might need some caution for widely divergent cell types / library sizes
def scale_vectors(density_vectors):
    overall_avg = np.mean([np.mean(density_vectors[k]) for k in density_vectors])
    mean_values = {d:np.mean(density_vectors[d]) for d in density_vectors}
    overall_avg = np.mean([mean_values[x] for x in mean_values])
    return {d:density_vectors[d]*overall_avg/mean_values[d] for d in mean_values}
    
def get_geomean(density_vectors):
    pseudocounts = {d:np.min(density_vectors[d][density_vectors[d]>0]) for d in density_vectors}
    
    gm_means = np.zeros(len(density_vectors[density_vectors.keys()[0]]))
    for k in density_vectors:
        gm_means += np.log(density_vectors[k]+pseudocounts[k])
    gm_means = gm_means / float(len(density_vectors))
    return np.exp(gm_means)
    
def get_mean(density_vectors):
    
    means = np.zeros(len(density_vectors[density_vectors.keys()[0]]))
    for k in density_vectors:
        means += density_vectors[k]
    means /= float(len(density_vectors))
    return means

def get_peak_subset(avg_values,numbers,density_vectors,correlation_limit):
    for i in range(int(np.max(numbers))+1):
        over = numbers>=i
        avg_cor = np.mean([spearmanr(avg_values[over],density_vectors[d][over])[0]
                    for d in density_vectors])
        if avg_cor > correlation_limit:
            break
    if i == np.max(numbers):
        print 'Check data: Individual Samples may be poorly captured by the mean'
    return np.where(numbers>=i)[0]

def get_peak_numbers(called_vectors):
    numbers = np.zeros(len(called_vectors[called_vectors.keys()[0]]))
    for c in called_vectors: numbers += called_vectors[c].astype(int)
    return numbers
    
def loess_geomean_norm(density_vectors, called_vectors, sample_number=SAMPLE_NUMBER,sample_method = 'random',correlation_limit = CORRELATION_LIMIT,
                        cv_number = 5,pre_normalize = True):

    numbers = get_peak_numbers(called_vectors)
    if pre_normalize:
        size_normed = scale_vectors(density_vectors)
    else:
        size_normed = density_vectors
    avg_values = get_mean(size_normed)
    
    pseudocounts = {d:np.min(size_normed[d][size_normed[d]>0]) for d in size_normed}
    gm_pseudo = np.mean([pseudocounts[x] for x in pseudocounts])
    
    decent_peaks = get_peak_subset(avg_values,numbers,size_normed,correlation_limit)

    sampled = select_uniform(avg_values[decent_peaks],number = sample_number,
                            sample_method = sample_method)
    sampled = decent_peaks[sampled]
    xvalues = np.log(avg_values+gm_pseudo)
    delta = np.percentile(avg_values,99)*DELTA_FRACTION
    difs = {k:np.log(size_normed[k]+pseudocounts[k])-np.log(avg_values+gm_pseudo) for 
                k in size_normed}
    ### don't need to cross validate so many times
    cv_set = SEED.choice(size_normed.keys(),size = min(cv_number,len(size_normed.keys())),replace = False)
    cv_fraction = np.mean([choose_fraction_cv(xvalues,difs[k],sampled,delta = delta)
                    for k in cv_set])
    
    normed_vectors = {}
    for k in size_normed:
        dif = difs[k]
        new_values = sm.lowess(dif[sampled],xvalues[sampled],
                                       return_sorted=False, it=4, frac = cv_fraction,
                                       delta = delta)
        interpolated2 = return_interpolated(xvalues, new_values, sampled)
        normed_vectors[k] =  size_normed[k]/np.exp(interpolated2)

    return normed_vectors

### creates size factors from the normalized counts for DESeq
def create_normalization_file(normed_vectors,raw_vectors,name_order,normalization_file = '',
                                header = False):
    size_matrix = []
    for i in name_order:
        n = normed_vectors[i]
        cp = raw_vectors[i]

        scales = cp/n
        potential_issue = np.logical_or(n<1,cp<1)
        scales[(n==0) | (cp == 0)] = 1 ## hack just in case some zero values, where size factor is not well defined
        size_matrix.append(scales)

    size_matrix = np.array(size_matrix)
    size_matrix[np.isnan(size_matrix)]=1.0
    size_matrix[np.logical_not(np.isfinite(size_matrix))]=1.0

    ### normalize size factors geometric mean
    geom_means = np.exp(np.mean(np.log(size_matrix)))
    size_matrix /= geom_means
    if normalization_file == '':
        n,out_file = tempfile.mkstemp()
        print 'Your Normalization File is %s' %out_file
    else:
        out_file = normalization_file
    if not header:
        np.savetxt(out_file,np.transpose(size_matrix),'%.5e',delimiter = '\t')
    else:
        np.savetxt(out_file,np.transpose(size_matrix),'%.5e',header = '\t'.join(name_order),comments = '',
                        delimiter = '\t')
    return out_file


def run_normalization_file(count_file, calls_file, geomean, output_file=None):
    density_data = np.loadtxt(count_file)
    calls_data = np.loadtxt(calls_file)
    calls_data = calls_data.astype(bool)
    norm_number = [1 / (np.sum(density_data[:, i]) * 2 / 100000.) for i in range(density_data.shape[1])]

    density_vectors = {i: density_data[:, i] * norm_number[i] for i in range(density_data.shape[1])}
    calls_vectors = {i: calls_data[:, i] for i in range(calls_data.shape[1])}

    if geomean:
        normed_vectors = loess_geomean_norm(density_vectors, calls_vectors,pre_normalize = False)
    else:
        normed_vectors = loess_group_norm(density_vectors, calls_vectors)
    base_name = os.path.split(count_file)[1].strip('.txt')
    base_dir = os.path.split(count_file)[0]

    orders = range(calls_data.shape[1])
    if output_file == None:
        output_file_factors = os.path.join(base_dir, base_name + '.size_factors')
        output_file_counts = os.path.join(base_dir, base_name + '.norm_counts')
    else:
        output_file_factors = output_file.rstrip('.txt') + '.size_factors'
        output_file_counts = output_file.rstrip('.txt') + '.norm_counts'

    density_vectors = {i: density_data[:, i] for i in orders}
    normed_vectors = {i: normed_vectors[i] / np.mean(norm_number) for i in orders}
    create_normalization_file(normed_vectors, density_vectors, orders,
                                  normalization_file=output_file_factors)
    np.savetxt(output_file_counts, np.transpose(np.array([normed_vectors[r] for r in orders])), fmt='%.6e')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process Count Files for Normalization.')
    parser.add_argument('counts_file', help='file of counts')
    parser.add_argument('calls_file', help='file of peak calls')
    parser.add_argument('--geomean', help='normalize to geometric mean',
                        default=True, action='store_true')
    parser.add_argument('-o', '--output', type=str, help='output file',
                        default=None)

    args = parser.parse_args()

    try:
        run_normalization_file(args.counts_file, args.calls_file, args.geomean, args.output)

    except Exception as e:
        log.error(e)