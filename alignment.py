import re
import os
import numpy as np
import pandas as pd
from itertools import compress
import math
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
# https://tiefenauer.github.io/blog/smith-waterman/

def threshold_and_get_rates(positives,negatives,subset_index):
    '''
    get the threshold and those that are then indicated as true and false positives
    '''

    threshold = sorted(positives)[subset_index]
    false_pos = get_rates(threshold,negatives)
    true_pos = get_rates(threshold,positives)
    return false_pos,true_pos

def get_score_dfs(gap_ext,gap,m,pos_pairs,neg_pairs,threshold):
    '''
    make a df from the scoring matrices
    '''
    pos_scores,pos_min,pos_align = run_sw_and_split_out(pos_pairs,m,gap,gap_ext)
    neg_scores,neg_min,neg_align = run_sw_and_split_out(neg_pairs,m,gap,gap_ext)

    # the index of the value from which to subset the sorted true positive list
    subset_index = int(len(pos_scores)-(len(pos_scores)*threshold/100)-1)

    # get the false and true positive rates, and normalized
    false_pos,true_pos = threshold_and_get_rates(pos_scores,neg_scores,subset_index)
    false_norm,true_norm = threshold_and_get_rates(pos_min,neg_min,subset_index)

    # get the gap, gap extention, false and true positives
    pd_stats = [true_pos,false_pos,gap,gap_ext,m,true_norm,false_norm]

    # collect the scores and other values for making roc curve
    all_score = pos_scores + neg_scores
    all_norm = pos_min + neg_min
    all_matrix = make_columns(m,m,pos_scores,neg_scores)
    all_labels = make_columns(1,0,pos_scores,neg_scores)
    all_gaps = make_columns(gap,gap,pos_scores,neg_scores)
    all_ext = make_columns(gap_ext,gap_ext,pos_scores,neg_scores)
    all_align = make_columns(pos_align,neg_align,pos_scores,neg_scores)
    pd_scores = pd.DataFrame({'scores':all_score,'norm_scores':all_norm,
        'matrix':all_matrix,'labels':all_labels,'gap':all_gaps,'ext':all_ext})#,'alignment':all_align

    return pd_stats,pd_scores

def run_sw_and_split_out(pairs,m,gap,gap_ext):
    '''
    run the smith waterman alg and split the scoring and source matrices into
    two outputs
    '''
    out = [smith_waterman(x.split()[0], x.split()[1],m,gap,gap_ext) for x in pairs]
    return [x[0] for x in out],[x[1] for x in out],[x[2] for x in out]

def get_rates(threshold,scores):
    '''
    get the false and true positive rates
    '''
    return sum([x > threshold for x in scores])/len(scores)

def get_norm_rates(threshold,min):
    '''
    get the normalized false and true positive rates
    '''
    return sum(list(compress(min, [x > threshold for x in min])))/len(min)

def make_columns(value_a,value_b,length_a,length_b):
    '''
    make columns for output data frame
    '''
    return [value_a]*len(length_a) + [value_b]*len(length_b)

def make_empty_matrices(length,width):
    '''
    make empty matrices populated with zeros
    '''
    return np.zeros((length,width))

def get_max(list):
    '''
    a list from which to get the max value
    '''
    return np.max(np.array(list))

def get_max_index(list):
    '''
    get the index of the max value in list
    '''
    return list.index(get_max(list))

def make_scoring_matrix(cost,str_a,str_b,open,extend):
    '''
    make a matrix of dimensions string_a x string_b
    populated by the match score between the two strings
    '''
    # get the lengths of the strings
    len_a = len(str_a) + 1
    len_b = len(str_b) + 1

    # initialze matrices with 0's
    # first is aligned, second is misaligned with gaps in a, third is misaligned with gaps in b
    M,A,B = make_empty_matrices(len_b, len_a),make_empty_matrices(len_b, len_a),make_empty_matrices(len_b, len_a)

    # initalize matrices for the sequence of commands
    # that indicate where each position came from
    M_c,A_c,B_c = make_empty_matrices(len_b, len_a),make_empty_matrices(len_b, len_a),make_empty_matrices(len_b, len_a)

    # # where there is no sequence to align to set to inf
    M[[i for i in range(1,len_b)],0] = -(math.inf)
    M[0,[j for j in range(1,len_a)]] = -(math.inf)
    B[[i for i in range(1,len_b)],0] = -(math.inf)
    A[0,[j for j in range(1,len_a)]] = -(math.inf)

    # gaps in first seq for every position in second for alignment to zero characters
    penalty = extend + open
    for j in range(1,len_a):
        B[0,j] = j * - penalty
    for i in range(1,len_b):
        A[i,0] = i * - penalty

    # for row and column along matrix
    for i, j in itertools.product(range(1, M.shape[0]), range(1, M.shape[1])):

        # the cost for A and B is going to depend on where the previous value is coming from
        # the cost for M is always going to be the cost matrix
        get_m = [M[i-1,j-1],A[i-1,j-1],B[i-1,j-1]]
        get_a = [M[i-1,j],A[i-1,j]-extend,B[i-1,j]]
        get_b = [M[i,j-1],A[i,j-1],B[i,j-1]-extend]
        a_max = get_max(get_a)
        b_max = get_max(get_b)
        m_max = get_max(get_m)
        A[i,j] = (0 if a_max < 0 else a_max) - open
        B[i,j] = (0 if b_max < 0 else b_max) - open
        M[i,j] = cost.loc[str_a[j-1],str_b[i-1]] + (0 if m_max < 0 else m_max)

        # get the index of the max value for each movement in the matrices
        # if the index is: 0 - source M, 1 - source A, 2 - source B
        A_c[i,j] = get_max_index(get_a)
        B_c[i,j] = get_max_index(get_b)
        M_c[i,j] = get_max_index(get_m)

    # get the matrix that has the highest value
    max_values = [M.max(),A.max(),B.max()]
    # get the index of the list housing the max elements from the matrices
    max_index = get_max_index(max_values)
    score_matrices = [M,A,B]
    source_matrices = [M_c,A_c,B_c]
    # get the score and source matrices associated that gives where that value came from
    score_matrix = score_matrices[max_index]
    source_matrix = source_matrices[max_index]

    # get the index of the max value in the score matrix
    max_instances = np.where(score_matrix == np.amax(score_matrix))
    pair_coords = list(zip(max_instances[0], max_instances[1]))

    # return only the last instance of the max value in the scoring matrix
    return score_matrix,source_matrix,pair_coords[-1]

def read_fasta(file):
    '''
    read in the fasta sequence and return sequence and protein name
    '''
    # get protein name
    name_a = file.split('.')[0]
    name_b = name_a.split('/')[1]

    # open file and read all lines, remove first for header
    with open(file,'r') as f:
        lines = f.read().splitlines()
        lines = lines[1:]

    # collect sequences
    sequence = []
    for line in lines:
        sequence = sequence + list(line)

    return sequence, name_b

def scoring_matrix_heatmap(matrix,str_a,str_b,name_a,name_b,start_position,align,score):
    '''
    take the scoring matrix and make heatmap
    '''
    sns.heatmap(matrix,xticklabels=False,yticklabels=False,vmin=0)
    plt.xlabel('String A, {0}'.format(name_a))
    plt.ylabel('String B, {0}'.format(name_b))
    plt.title('Starting:{0}, score:{1}'.format(start_position,score))
    plt.suptitle('{0}'.format(align))

    # save plot to image directory
    outdir = pathlib.Path('images')
    outfile = outdir / 'Heatmap_{0}_{1}.png'.format(name_a,name_b)
    outdir.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(outfile),format='png')
    plt.clf()
    plt.close()

def optimal_traceback(score_matrix, string_a, string_a_='',previous_i=0):
    '''
    recurisively for each substring in string a
    find the optimal traceback path within the scoring matrix
    get each last max occurance starting from the bottom of the scoring matrix

    doesn't yet account for direction recieved!!
    '''
    # invert the matrix in order to get the last value from both axes
    score_matrix_flip = np.flip(np.flip(score_matrix, 0), 1)

    # get the index of the max value in the score matrix
    i_, j_ = np.unravel_index(score_matrix_flip.argmax(), score_matrix_flip.shape)

    # adjust the index values
    i, j = np.subtract(score_matrix.shape, (i_ + 1, j_ + 1))

    # if 0 it means you've gone through the entire matrix
    if score_matrix[i, j] == 0:
        return string_a_, j

    # other wise add to the string value or add gap string
    string_a_ = string_a[j - 1] + '-' + string_a_ if previous_i - i > 1 else string_a[j - 1] + string_a_

    # recursivley perform through axis b string
    return optimal_traceback(score_matrix[0:i, 0:j], string_a, string_a_, i)

def read_in_cost_matrix(matrix_file):
    '''
    read in the costmatrix
    '''
    # get the rows that start with comment in the cost matrix
    exclude = [i for i, line in enumerate(open(os.path.join('matrices',matrix_file),"r")) if line.startswith('#')]

    # open matrix file
    m_file = open(os.path.join('matrices', matrix_file), "r")

    # read in the matrix, skipping the rows that are comments
    cost_matrix =  pd.read_fwf(m_file,skiprows = exclude[0:])

    # set the index to be the same as the column names
    cost_matrix.set_index(cost_matrix.columns.values,inplace=True)
    return cost_matrix

def smith_waterman(file_a,file_b,matrix_file,open_gap,extend_gap):
    '''
    perform smith waterman algorithm, getting the traceback along the
    matrix of scored values between strings a and b
    '''

    # get the sequence for files a and b
    string_a, name_a = read_fasta(file_a)
    string_b, name_b = read_fasta(file_b)

    # read in cost matrix
    cost_matrix = read_in_cost_matrix(matrix_file)

    # make sure strings are upper case
    string_a = [x.upper() for x in string_a]
    string_b = [x.upper() for x in string_b]

    # get the length of the smallest string pair
    min_length = min(len(string_a),len(string_b))

    # create the scoring matrix for strings a and b
    score_matrix,source_matrix,start_position = make_scoring_matrix(cost_matrix,string_a,string_b,open_gap,extend_gap)

    # get the optimal trace back through the scoring matrix and find the starting point in the string
    alignment, start = optimal_traceback(score_matrix, string_a)
    # print(alignment)

    scoring_matrix_heatmap(score_matrix,string_a,string_b,name_a,name_b,start_position,alignment,score_matrix[start_position])

    # get the score of the matrix
    return score_matrix[start_position],score_matrix[start_position]/min_length,alignment
