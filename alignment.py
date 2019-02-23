import re
import os
import numpy as np
import pandas as pd
import math
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
# import numba as nb
# https://tiefenauer.github.io/blog/smith-waterman/

# @nb.jit(parallel=True)
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
    M = np.zeros((len_b, len_a))
    A = np.zeros((len_b, len_a))
    B = np.zeros((len_b, len_a))

    # initalize matrices for the sequence of commands
    # that indicate where each position came from
    M_c = np.zeros((len_b, len_a))
    A_c = np.zeros((len_b, len_a))
    B_c = np.zeros((len_b, len_a))

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
        a_max = np.max(np.array(get_a))
        b_max = np.max(np.array(get_b))
        m_max = np.max(np.array(get_m))
        A[i,j] = (0 if a_max < 0 else a_max) - open
        B[i,j] = (0 if b_max < 0 else b_max) - open
        M[i,j] = cost.loc[str_a[j-1],str_b[i-1]] + (0 if m_max < 0 else m_max)

        # get the index of the max value for each movement in the matrices
        # if the index is: 0 - source M, 1 - source A, 2 - source B
        A_c = get_a.index(np.max(np.array(get_a)))
        B_c = get_b.index(np.max(np.array(get_b)))
        M_c = get_m.index(np.max(np.array(get_m)))

    # get the matrix that has the highest value last element
    max_values = [M[-1,-1],A[-1,-1],B[-1,-1]]
    max_index = max_values.index(np.max(np.array(max_values)))
    matrices = [M,A,B]
    source_matrices = [M_c,A_c,B_c]
    max = matrices[max_index]
    source = source_matrices[max_index]

    return max,source

def optimal_traceback(score_matrix, string_b, string_b_='',previous_i=0):
    '''
    recurisively for each substring in string b
    find the optimal traceback path within the scoring matrix
    get each last max occurance starting from the bottom of the scoring matrix
    '''
    # invert the matrix in order to get the last value from both axes
    score_matrix_flip = np.flip(np.flip(score_matrix, 0), 1)

    # get the index of the max value in the score matrix
    i_, j_ = np.unravel_index(score_matrix_flip.argmax(), score_matrix_flip.shape)

    # adjust the index values
    i, j = np.subtract(score_matrix.shape, (i_ + 1, j_ + 1))

    # if 0 it means you've gone through the entire matrix
    if score_matrix[i, j] == 0:
        return string_b_, j

    # other wise add to the string value or add gap string
    string_b_ = string_b[j - 1] + '-' + string_b_ if previous_i - i > 1 else string_b[j - 1] + string_b_
    # print(string_b, string_b_,previous_i)
    # print(i_, j_)
    # print(i,j)

    # recursivley perform through axis b string
    return optimal_traceback(score_matrix[0:i, 0:j], string_b, string_b_, i)

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

def scoring_matrix_heatmap(matrix,str_a,str_b,name_a,name_b):
    '''
    take the scoring matrix and make heatmap
    '''
    sns.heatmap(matrix,xticklabels=False,yticklabels=False,vmin=0)
    plt.xlabel('String A, {0}'.format(name_a))
    plt.ylabel('String B, {0}'.format(name_b))
    # plt.title('Heatmap for {0} and {1}'.format(name_a,name_b))

    # save plot to image directory
    outdir = pathlib.Path('images')
    outfile = outdir / 'Heatmap_{0}_{1}.png'.format(name_a,name_b)
    outdir.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(outfile),format='png')
    plt.clf()
    plt.close()

def smith_waterman(file_a,file_b,matrix_file,open_gap,extend_gap):
    '''
    perform smith waterman algorithm, getting the traceback along the
    matrix of scored values between strings a and b
    '''

    # get the sequence for files a and b
    string_a, name_a = read_fasta(file_a)
    string_b, name_b = read_fasta(file_b)

    # get the rows that start with comment in the cost matrix
    exclude = [i for i, line in enumerate(open(os.path.join('matrices',matrix_file),"r")) if line.startswith('#')]

    # open matrix file
    m_file = open(os.path.join('matrices', matrix_file), "r")

    # read in the matrix, skipping the rows that are comments
    cost_matrix =  pd.read_fwf(m_file,skiprows = exclude[0:])

    # set the index to be the same as the column names
    cost_matrix.set_index(cost_matrix.columns.values,inplace=True)

    # make sure strings are upper case
    string_a = [x.upper() for x in string_a]
    string_b = [x.upper() for x in string_b]

    # get the length of the smallest string pair
    min_length = min(len(string_a),len(string_b))

    # create the scoring matrix for strings a and b
    score_matrix,source_matrix = make_scoring_matrix(cost_matrix,string_a,string_b,open_gap,extend_gap)

    scoring_matrix_heatmap(score_matrix,string_a,string_b,name_a,name_b)

    # get the optimal trace back through the scoring matrix and find the starting point in the string
    # string_b_, position = optimal_traceback(score_matrix, string_b)

    # return the starting point in the string, and the length of the matched string b
    #position, position + len(string_b_)

    # get the score of the matrix
    return score_matrix[-1,-1],score_matrix[-1,-1]/min_length
