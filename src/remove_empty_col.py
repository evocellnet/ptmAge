#!/usr/bin/env python2

import sys
from Bio import SeqIO


# Load sequences
matrix = {}
matrix_id = []
input_handle = open(sys.argv[1], "r")
for record in SeqIO.parse(input_handle, "fasta"):
    matrix[record.id] = record.seq
    matrix_id.append(record.id)
    
columns = []
keys = matrix.keys()
for i in matrix[keys[0]]:
    columns.append(0)
#print len(columns)

#print matrix


# Keep "X"
for i in matrix:
    for j in range(len(matrix[i])):
        #print str(j)+" "+matrix[i][j]
        aa = matrix[i][j]
        if aa == "-":
            columns[j] = columns[j]+0
        #elif aa == "X":
        #    columns[j] = columns[j]+0     
        else:
            columns[j] = columns[j]+1

matrix_new = {}

tag = 0
for i in matrix:
    tmp = []
    j = 0
    for pos in range(len(columns)):
        if columns[pos] != 0:
            j = j+1
          ##   if tag == 0:
##                 print j,"\t",pos+1
            #print matrix[i][pos]
            tmp.append(matrix[i][pos])
        #else:
            #if tag == 0:
                #print pos
    matrix_new[i] = tmp
    tag = 1
#print matrix2

# Write file
file_out = open(sys.argv[2],"w")
for gene in matrix_id:
    file_out.write(">"+gene+"\n")
    file_out.write("".join(matrix_new[gene])+"\n")
file_out.close()
