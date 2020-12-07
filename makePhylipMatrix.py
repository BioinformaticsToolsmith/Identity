#
# Author: Hani Z. Girgis, PhD
#
#
# Purpose: This program converts all-versus-all identity scores produced 
# by Identity to Phylip distance matrix
# Requirements: Python 3 with NumPy.
#

import numpy as np
import sys

   
def make_matrix(file_in, file_out):    
    file = open(file_in, 'r')
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
    
    # Identity does not report the sequence versus itself
    for i in range(0, label_count):
        matrix[i][i] = 0.00000000
            
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
        matrix[id_1][id_2] = distance
        matrix[id_2][id_1] = distance
      
    file.close()
 
    # Write Phylip matrix
    sorted_key_list = sorted(label_dict.keys(), key=lambda key: label_dict[key])
    with open(file_out, 'w') as file_object:
        file_object.write(str(label_count) + '\n')
        for key in sorted_key_list:
            t = list()
            i = label_dict[key]
            for j in range(0, label_count):
                t.append(str.format('{0:.8f}', matrix[i][j])) #round(matrix[i][j], 8))
            s = ' '    
            t = [str(x) for x in t]
            # End a species name with \t
            file_object.write(key[1:] + '\t' + s.join(t) + '\n')

if len(sys.argv) != 3:
    print("Use: python3 ", sys.argv[0], "allVsAllIdentityFile outputFileName")
    print()
    print("Please provide an all-versus-all file produced by Identity and an output file name.")
    print()
    print("Example: python3 " + sys.argv[0] + " all_vs_all_identity_scores matrix.phylip")
    print()
else:
    make_matrix(sys.argv[1], sys.argv[2])
