#!/usr/bin/python

import os
import sys

# Define version
version = sys.argv[1].split(".")[-4]

#print sys.argv[1].split(".")
#print version

#sys.exit()

tree_directory = sys.argv[2]

# Load mega tree file Compara.XX.protein.aa.fasta
alignment_file = open(sys.argv[1], "r")

#gene_to_aa_file = open("Gene_to_aa_file.txt", "w")

size_name = -1
# if len(sys.argv) > 2:
#     size_name = int(sys.argv[3])


tree_number = 1
tree_name = "ENSTREE_"+str(tree_number).zfill(5)
supra_folder = "ENSTREE_"+str(int(str(tree_number).zfill(5)[0:2])+1).zfill(2)+"000"
folder = "ENSTREE_"+str(tree_number).zfill(5)

os.system("mkdir -p "+tree_directory+"/"+supra_folder+"/"+folder)

aa_file = open(tree_directory+"/"+supra_folder+"/"+folder+"/"+tree_name+".aa.fasta", "w")

i = 0
while 1:
    line = alignment_file.readline()
    i = i+1
    if line == "":
        break
    if line[0:2] == "//":
        aa_file.close()                                  # We close the previous file.
        tree_number = tree_number+1                      # We increase by +1.
        tree_name = "ENSTREE_"+str(tree_number).zfill(5) #
        supra_folder = "ENSTREE_"+str(int(str(tree_number).zfill(5)[0:2])+1).zfill(2)+"000"
        folder = "ENSTREE_"+str(tree_number).zfill(5)
        os.system("mkdir -p "+tree_directory+"/"+supra_folder+"/"+folder)
        aa_file = open(tree_directory + "/"+supra_folder+"/"+ \
                       folder+"/"+tree_name+".aa.fasta", "w") # We open a new file to write in.
        #print tree_number, "ENSTREE_"+str(tree_number).zfill(5)
        #sys.exit()
    elif len(line) >= 1:
        aa_file.write(line)

aa_file.close()
alignment_file.close()
