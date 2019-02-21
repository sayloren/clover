import os
import pandas as pd
import re
from .alginment import smith_waterman

# 1) get pairs
pos_pairs = [line.strip() for line in 'pairs/Pospairs.txt']
neg_pairs = [line.strip() for line in 'pairs/Nospairs.txt']

# initialze the collection data frame
collect = []

# 2) get matrix scoring fasta pair for each matrix
# for m in os.listdir('matrices'):
# for gap in range(1,21):
# for ext in range(1,6):
# pos_scores = [smith_waterman(x[0], x[1],m,gap,gap_ext) for x in pos_pairs]
# neg_scores = [smith_waterman(x[0], x[1],m,gap,gap_ext) for x in neg_pairs]
pos_scores = [smith_waterman(x[0], x[1],'matrices/BLOSUM50',3,2) for x in pos_pairs]
neg_scores = [smith_waterman(x[0], x[1],'matrices/BLOSUM50',3,2) for x in neg_pairs]
print(pos_scores)

#         #Threshold for 0.7 True Positive Rate. (there are 50 Pos_pairs)
#         threshold = sorted(pos_scores)[14]
#         false_pos = np.sum(np.array(neg_scores) > threshold)/len(neg_scores)
#         true_pos = np.sum(np.array(pos_scores) > threshold)/len(pos_scores)
#         fp_result.append([gap,ext,false_pos,true_pos])
