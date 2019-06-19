import argparse
import logging
import numpy as np
import os
import time
from scipy.stats import expon

import bedtools
import smoothers_lowess as sm

# script log created as deseq2_{today}.log in working directory
today = time.strftime("%d%m%Y")
FORMAT = "[%(name)s | Line:%(lineno)s| Function:%(funcName)s ] %(message)s"
logging.basicConfig(filename="deseq2_{}.log".format(today), level=logging.INFO, filemode='w',
                    format=FORMAT)
log = logging.getLogger('normalization_script')


SEED_NUMBER = 1832245
SEED_ITERATION = 10000
SEED = np.random.RandomState(SEED_NUMBER)


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


def update_seed():
    global SEED_NUMBER
    global SEED

    direction = 1
    if SEED_NUMBER > 4294000000:
        direction = -1
    elif SEED_NUMBER < 100000:
        direction = 1
    SEED_NUMBER += direction * (SEED_ITERATION + SEED.randint(1, SEED_ITERATION))
    SEED.seed(SEED_NUMBER)

def filter_on_peaks(core_subset,called1,called2):###
    good = np.intersect1d(np.where(called1>0)[0],np.where(called2>0)[0])
    return np.intersect1d(core_subset,good)

def return_ordered_indexes(density_vectors,called_vectors):
    counts = np.zeros(len(density_vectors[density_vectors.keys()[0]]))
    gm_means = np.ones(len(counts))
    for k in density_vectors:
        gm_means[called_vectors[k]] *= density_vectors[k][called_vectors[k]]
        counts[called_vectors[k]] += 1
    gm_means = gm_means ** (1/counts)
    sorted_indexes = np.argsort(gm_means)
    return sorted_indexes,gm_means,counts

def filter_by_range(sorted_list,threshold):
    LIMIT = 40000
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


def sample_from_range(sorted_geometric_means, sorted_indexes, percentile_cutoffs=[.05, 1],
                      sample_number=10000, tolerance_percentage=0.01, subset=None, top_cutoff=0.05, hard_cutoffs=None):
    cut_indexes = [int((len(sorted_geometric_means) - 1) * r) for r in percentile_cutoffs]
    while cut_indexes[-1] >= (len(sorted_geometric_means) - 1):
        cut_indexes[-1] -= 1
    cut_indexes[-1] = filter_by_range(sorted_geometric_means[0:cut_indexes[-1]], top_cutoff)

    if hard_cutoffs != None:
        cut_indexes[0] = np.searchsorted(sorted_geometric_means, hard_cutoffs[0])
        cut_indexes[1] = np.searchsorted(sorted_geometric_means, hard_cutoffs[1])
    sorted_log_heights = sorted_geometric_means
    cutoffs = [sorted_log_heights[i] for i in cut_indexes]

    tolerance = (cutoffs[1] - cutoffs[0]) * tolerance_percentage

    sample_slices = []
    indexes = set([])
    min_mean = sorted_log_heights[0]
    max_mean = sorted_log_heights[-1]
    for r in xrange(sample_number):
        hcenter = SEED.random_sample() * (cutoffs[1] - cutoffs[0]) + cutoffs[0]
        update_seed()
        sample_slices.append((max(min_mean, hcenter - tolerance), min(max_mean, hcenter + tolerance)))
    tolerable_min_index = np.searchsorted(sorted_log_heights, [t[0] for t in sample_slices], side='left')
    tolerable_max_index = np.searchsorted(sorted_log_heights, [t[1] for t in sample_slices], side='right') - 1

    for i, x in enumerate(tolerable_min_index):
        if x > tolerable_max_index[i]:
            continue
        if subset != None:
            acceptable = np.intersect1d(np.arange(x, tolerable_max_index[i] + 1), subset)
            if len(acceptable) == 0:
                continue
            indexes.add(int(SEED.choice(acceptable)))
            update_seed()
        else:
            acceptable = np.setdiff1d(np.arange(x, tolerable_max_index[i] + 1), indexes)
            if len(acceptable) == 0:
                continue
            indexes.add(int(SEED.choice(acceptable)))
            update_seed()
    index_random = sorted_indexes[list(indexes)]
    return index_random

def group_norm(density_vectors, called_vectors, method=None, data_type='DNAse',
               sample_number=10000, sample_once=False):
    sorted_indexes, gm_means, numbers = return_ordered_indexes(density_vectors, called_vectors)
    called_twice = np.where(numbers > 1)[0]

    if sample_once:
        sampled = sample_from_range(gm_means[sorted_indexes[called_twice]], sorted_indexes[called_twice],
                                    sample_number=sample_number,
                                    percentile_cutoffs=PEAK_THRESHOLDS[data_type], top_cutoff=TOP_THRESHOLDS[data_type])
    else:
        sampled = None

    pseudocount = 1.

    for k in density_vectors:
        pseudocount = min(pseudocount, np.min(density_vectors[k][np.where(density_vectors[k] > 0)[0]]))

    ref_distances = {}
    for r in density_vectors:
        ref_distances[r] = np.mean((density_vectors[r][called_vectors[r]] - gm_means[called_vectors[r]]) ** 2)
    ref_exp = sorted([r for r in ref_distances], key=lambda k: ref_distances[r])[0]

    normed_vectors = {}
    for k in density_vectors:
        if k == ref_exp:
            normed_vectors[k] = density_vectors[k]
        else:
            if FILTER_CALLED[data_type] and sampled != None:
                current_sample = filter_on_peaks(sampled, called_vectors[ref_exp], called_vectors[k])
            else:
                current_sample = sampled
            normed_vectors[k] = method(density_vectors[ref_exp], density_vectors[k], called_vectors[ref_exp],
                                       called_vectors[k], sampled=current_sample,
                                       filter_on_logratio=TMM_RATIOS[data_type], pseudocount=pseudocount)
    return normed_vectors

def extrapolate(interpolated,to_extrap,base,predict):
    EXTRAPOLATION_FRACTION = 0.1
    sample_order = np.argsort(base)
    min_index_range = (sample_order[0],sample_order[int(EXTRAPOLATION_FRACTION*len(sample_order))])
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

def lowess_norm_pair(peaks1, peaks2, called1, called2, sampled=None, pseudocount=1):
    geomean = np.sqrt(peaks1 * peaks2)
    if sampled == None:
        called_both = np.where(called1 * called2 > 0)[0]
        peaks1_s = np.argsort(peaks1[called_both])
        peaks2_s = np.argsort(peaks2[called_both])
        sampled_1 = sample_from_range(np.sqrt(peaks1[called_both][peaks1_s]),
                                      peaks1_s, percentile_cutoffs=[0, 1], sample_number=15000, top_cutoff=0.025)
        sampled_2 = sample_from_range(np.sqrt(peaks2[called_both][peaks2_s]),
                                      peaks2_s, percentile_cutoffs=[0, 1], sample_number=15000, top_cutoff=0.025)

    uncalled_fraction = np.where(called1 * called2 == 0)[0]

    min_1 = np.min(peaks1[called_both][sampled_1])
    max_1 = np.max(peaks1[called_both][sampled_1])
    number_1 = int(len(sampled_1) * min_1 / float(max_1 - min_1))
    min_2 = np.min(peaks2[called_both][sampled_2])
    max_2 = np.max(peaks2[called_both][sampled_2])
    number_2 = int(len(sampled_2) * min_2 / float(max_2 - min_2))

    if len(uncalled_fraction) <= number_1:
        uncalled_sample1 = np.arange(len(uncalled_fraction))
        print 'Warning: Few Uncalled Peaks'
    else:
        uncalled_sample1 = SEED.choice(len(uncalled_fraction), number_1, replace=False)
        update_seed()
    if len(uncalled_fraction) <= number_2:
        print 'Warning: Few Uncalled Peaks'
        uncalled_sample2 = np.arange(len(uncalled_fraction))
    else:
        uncalled_sample2 = SEED.choice(len(uncalled_fraction), number_2, replace=False)
        update_seed()
    #### need to add in all the way down to zero
    sampled_1 = np.append(called_both[sampled_1], uncalled_fraction[uncalled_sample1])
    sampled_2 = np.append(called_both[sampled_2], uncalled_fraction[uncalled_sample2])


    new1 = sm.lowess(geomean[sampled_1], peaks1[sampled_1],
                                   return_sorted=False, it=4, frac=0.4, delta=0.01 * max(peaks1[sampled_1]))
    new2 = sm.lowess(geomean[sampled_2], peaks2[sampled_2],
                                   return_sorted=False, it=4, frac=0.4, delta=0.01 * max(peaks2[sampled_2]))

    interpolated2 = return_interpolated(peaks2, new2, sampled_2, pseudocount)

    sampled_order = np.argsort(new1)
    interpolated3 = np.interp(interpolated2, new1[sampled_order], peaks1[sampled_1][sampled_order])
    interpolated3 = extrapolate(interpolated3, interpolated2, new1, peaks1[sampled_1])

    return interpolated3


def loess_group_norm(density_vectors,called_vectors,data_type = 'DNAse'):
    normed_values = group_norm(density_vectors,called_vectors,method=lowess_norm_pair,
                    data_type=data_type,sample_number=25000,sample_once = False)
    return normed_values


def loess_geomean_norm(density_vectors, called_vectors, sample_number=30000, data_type='DNAse'):
    #### because to the way this is set up, it only works if a
    #### the same general peaks are the same:  look at called twice
    ### might want to change it at some point
    #### also if density vectors are not normalized first, run into problems

    numbers = np.zeros(len(density_vectors[density_vectors.keys()[0]]))
    gm_means = np.ones(len(numbers))
    pseudocount = np.min([np.min(density_vectors[d]) for d in density_vectors])
    for k in density_vectors:
        gm_means *= (density_vectors[k] + pseudocount)
        numbers[called_vectors[k]] += 1
    gm_means = gm_means ** (1 / float(len(density_vectors)))
    # sorted_indexes = np.argsort(gm_means)
    called_twice = np.where(numbers > 1)[0]

    normed_vectors = {}
    for k in density_vectors:
        vector_sort = np.argsort(density_vectors[k][called_twice])
        sampled = sample_from_range(np.sqrt(density_vectors[k][called_twice][vector_sort]), called_twice[vector_sort],
                                    sample_number=sample_number,
                                    percentile_cutoffs=PEAK_THRESHOLDS[data_type], top_cutoff=TOP_THRESHOLDS[data_type])
        sampled = sampled[gm_means[sampled] != 0]
        uncalled_fraction = np.where(called_vectors[k] == 0)[0]

        min_gm = np.sqrt(np.min(density_vectors[k][sampled]))
        max_gm = np.sqrt(np.max(density_vectors[k][sampled]))
        number_gm = int(len(sampled) * min_gm / float(max_gm - min_gm))

        if len(uncalled_fraction) <= number_gm:
            uncalled_sample_gm = np.arange(len(uncalled_fraction))
            print('Warning: Few Uncalled Peaks')
        else:
            uncalled_sample_gm = SEED.choice(len(uncalled_fraction), number_gm, replace=False)
            update_seed()
        sample_subset = np.append(sampled, uncalled_fraction[uncalled_sample_gm])

        pseudocount = min(pseudocount, np.min(density_vectors[k][np.where(density_vectors[k] > 0)[0]]))
        new1 = sm.lowess(gm_means[sample_subset], density_vectors[k][sample_subset],
                                       return_sorted=False, it=4, frac=0.4,
                                       delta=0.01 * max(density_vectors[k][sample_subset]))

        interpolated = return_interpolated(density_vectors[k], new1, sample_subset, pseudocount)
        normed_vectors[k] = interpolated

    return normed_vectors

def create_normalization_file(normed_vectors,raw_vectors,name_order,normalization_file = '',
                                header = False):
    size_matrix = []
    for i in name_order:
        n = normed_vectors[i]
        cp = raw_vectors[i]
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


def run_normalization_file(count_file, calls_file, geomean, output_file=None):
    density_data = np.loadtxt(count_file)
    calls_data = np.loadtxt(calls_file)
    calls_data = calls_data.astype(bool)
    norm_number = [1 / (np.sum(density_data[:, i]) * 2 / 100000.) for i in range(density_data.shape[1])]

    density_vectors = {i: density_data[:, i] * norm_number[i] for i in range(density_data.shape[1])}
    calls_vectors = {i: calls_data[:, i] for i in range(calls_data.shape[1])}

    if geomean:
        normed_vectors = loess_geomean_norm(density_vectors, calls_vectors)
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
                        default=False, action='store_true')
    parser.add_argument('-o', '--output', type=str, help='output file',
                        default=None)

    args = parser.parse_args()

    try:
        run_normalization_file(args.counts_file, args.calls_file, args.geomean, args.output)

    except Exception as e:
        log.error(e)