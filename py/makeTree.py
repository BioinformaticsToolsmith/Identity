#
# Author: Hani Z. Girgis, PhD
# Purpose: This program performs hierarchical clustering using the 
#    all-versus-all identity scores produced by Identity.
# Requirements: Python 3 with SciPy and Matplotlib packages.
#

import numpy as np
import sys

from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import leaves_list
import matplotlib.pyplot as plt

   
def generate_tree(file, displayTree):    
    file = open(file, 'r')
    # List of labels
    label_dict = {}
    seq_index = 0
    for line in file:
        token_list = line.strip().split("\t")
         
        if not token_list[0] in label_dict:
            label_dict[token_list[0]] = seq_index
            seq_index += 1
         
        if not token_list[1] in label_dict:    
            label_dict[token_list[1]] = seq_index
            seq_index += 1
 
    label_count = len(label_dict.keys())
 
    # The matrix 
    matrix = np.ones((label_count, label_count))
    # Go to the begining of the file
    file.seek(0)
    for line in file:
        token_list = line.strip().split("\t")
        distance = 1.0 - float(token_list[2])
 
        if(distance < 0.0):
            distance = 0.0
        elif(distance > 1.0):
            distance = 1.0
 
        id_1 = label_dict[token_list[0]]
        id_2 = label_dict[token_list[1]]
        if(id_1 < id_2):
            matrix[id_1][id_2] = distance
        else:
            matrix[id_2][id_1] = distance
      
    file.close()
 
    key_list = label_dict.keys()
    name_list = [-1] * len(key_list)
    for k in key_list:
        name_list[label_dict[k]] = k; 
     
    l = list(matrix[np.triu_indices(label_count, 1)])
    row_clusters = linkage(l, method='average', optimal_ordering=True)
     
    # Print ordered leaves 
    if displayTree == 0:
        leaf_list = leaves_list(row_clusters)
        for index in leaf_list:
            print(name_list[index])
    else:
        plt.figure()
        dendrogram(row_clusters, orientation='left', labels=name_list, leaf_font_size=5, distance_sort='descending', show_leaf_counts=True)
        plt.tight_layout()
        plt.show()


if len(sys.argv) != 3:
    print("Use:", sys.argv[0], "allVsAllIdentityFile sortedLeavesOrTree")
    print()
    print("Please provide an all-versus-all file produced by Identity and 0 (display sorted leaves) or 1 (display a tree).")
    print()
    print("Display only sorted leaves: python3 makeTree.py all_vs_all_identity_scores 0")
    print("Display tree: python3 makeTree.py all_vs_all_identity_scores 1")
    print()
else:
    generate_tree(sys.argv[1], int(sys.argv[2]))
