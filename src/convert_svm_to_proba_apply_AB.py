#!/usr/bin/env python3

import sys
import math

#def sigmoid(x):
#  return 1 / (1 + math.exp(-x))


### Input parameters

# out = array of SVM outputs
# target = array of booleans: is ith example a positive example
# prior1 = number of positive examples
# prior2 = number of negative examples

### Output parameteres
# A, B = parameters of sigmoid


# Define variable


prior0 = 0.0
prior1 = 0.0


# Load SVM output
svm_output = []
svm_output_file = open(sys.argv[1], "r")
while 1:
    line = svm_output_file.readline()
    if line == "":
        break
    line = line.rstrip()
    svm_score = float(line)
    svm_output.append(svm_score)
svm_output_file.close()


# Load SVM target
target = []
gene_site_list = []
svm_target_file = open(sys.argv[2], "r")
while 1:
    line = svm_target_file.readline()
    if line == "":
        break
    tab = line.split()
    if tab[0] != "#":
        svm_target = tab[0]
        if svm_target[0] == "+":
            prior0 = prior0+1
            target.append(1)
        elif svm_target[0] == "-":
            prior1 = prior1+1
            target.append(0)
        #print tab
        gene_site = tab[-1][1:]
        #print gene_site
        gene_site_list.append(gene_site)
svm_target_file.close()

#sys.exit()

A, B = 0.0, 0.0

file_in = open(sys.argv[4], "r")
while 1:
    line = file_in.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split()
    A, B = float(tab[0]), float(tab[1])
file_in.close()

#print A, B
#sys.exit()


# Problem with Value domain math error
# => Use Sigmoid instead



# sys.exit()
output = open(sys.argv[3], "w")
i = 0
for score in svm_output:
    proba = 1/(1+math.exp(A*score+B))
    #proba = sigmoid(score)
    sign = "NEG"
    if target[i] == 1:
        sign = "POS"
    gene_site = gene_site_list[i]
    gene = "_".join(gene_site.split("_")[0:-1])
    site = gene_site.split("_")[-1]
    i = i+1
    # Write line
    output.write(gene+"\t"+site+"\t"+sign+"\t"+str(score)+"\t"+str(proba)+"\n")
output.close()
