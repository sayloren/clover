# clover  

algorithms assignment 3, sequence alignment  

[![Build
Status](https://travis-ci.org/sayloren/clover.svg?branch=master)](https://travis-ci.org/sayloren/clover)  

[Travis Build Results](https://travis-ci.org/sayloren/clover)  

:see_no_evil: :hear_no_evil: :speak_no_evil:  

## UCSF BMI 203 Algorithms Homework 2018  

#### Smith-Waterman Sequence Alignment Feb 22  

##### To run sort  
```
python __main__.py -t -g -e -s -i
```

-t the percentage at which to threshold true positives  
-g the gap opening size to run as a range from 1  
-e the gap extention penalty to run as a range from 1  
-s if want to run for a single set of gap and ext conditions rather than a range  
-i if data file should be read in, other wise the algorithm is run and the matrix generated  

##### To run sort tests  
```
python -m pytest  
```

##### The ROC for the best performing matrix  
![a](/images/ROC_Best_.png)  

##### The ROC normed for the best performing matrix  
![a](/images/ROC_Norm_.png)  

##### Example of heatmap for score matrix between two sequences  
![a](/images/Heatmap_prot-0004_prot-0008.png)  
