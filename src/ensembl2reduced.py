#!/usr/bin/env python

import os
import sys
from Bio import SeqIO

# Load gene families
families_list = []
families_file = open(sys.argv[1], "r")
while 1:
    line = families_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t")
    family = tab[2]
    families_list.append(family)
families_file.close()

# Make them unique
families_list = list(set(families_list))
families_list.sort()

#print families_list

# Define species to keep
species_list = []
species_file = open(sys.argv[2], "r")
while 1:
    line = species_file.readline()
    if line == "":
        break
    species = line.rstrip()
    species_list.append(species.replace(" ", "_"))
species_file.close()

#print species_list
#sys.exit()

species_tag = sys.argv[3]

iterator = 0
for path in families_list:
    iterator = iterator+1

    subpath = path+"/"+species_tag
    os.system("mkdir -p "+subpath)
    os.chdir(subpath)

    family = path.split("/")[-1]
    aln_file = "../"+family+".aa.fasta"
    tree_file = "../"+family+".tree"
    gene_file = "../"+family+".list"

    # Define subest of genes
    reduced_gene_set = []
    fine_in = open(gene_file, "r")
    while 1:
        line = fine_in.readline()
        if line == "":
            break
        line = line.rstrip()
        tab = line.split("\t")
        gene = tab[0]
        species = tab[1].replace(" ", "_")
        if species in species_list:
            reduced_gene_set.append(gene)
    fine_in.close()
    # print reduced_gene_set

    # Reduce alignment
    reduced_aln_file = family+"_reduced.aa.fasta"
    file_out = open(reduced_aln_file, "w")
    handle = open(aln_file, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        if record.id in reduced_gene_set:
            file_out.write(">"+record.id+"\n")
            file_out.write(str(record.seq)+"\n")
    handle.close()
    file_out.close()

    # Remove empty column
    script_file = open(family+".sh", "w")
    script_file.write("../../../../../src/remove_empty_col.py "+reduced_aln_file+" "+ \
                      family+"_reduced.aa.nogap.fasta\n")

    # Convert to phylip
    script_file.write("../../../../../src/convert_fasta2phylip.py "+family+"_reduced.aa.nogap.fasta"+ \
                      " "+family+"_reduced.aa.nogap.phy\n")

    # Write reduced .tree file
    script_file.write("nw_prune -v "+tree_file+" "+" ".join(reduced_gene_set)+ \
                      " > "+family+"_reduced.tree\n")

    # Recompute branch length
    script_file.write("phyml --quiet -i "+family+"_reduced.aa.nogap.phy -d aa -m LG -c 4 -a e -u " \
                      +family+"_reduced.tree -o lr -b 0\n")

    # Convert to ultrametric
    r_script = open(family+".ultra.r", "w")
    r_script.write("library(ape)\n")
    r_script.write("tree <- read.tree(\""+family+"_reduced.aa.nogap.phy_phyml_tree.txt\")\n")
    r_script.write("tree <- chronos(tree)\n")
    r_script.write("write.tree(tree, \""+family+"_reduced.aa.nogap.ultrametric.tree\")\n")
    r_script.close()

    script_file.write("R CMD BATCH "+family+".ultra.r\n")
    script_file.close()

    # Execute script
    os.system("bsub -M 10000 -R 'rusage[mem=10000]' -e /dev/null -o /dev/null sh "+family+".sh")

    os.chdir("../../../../")

    #if iterator == 5:
    #    sys.exit()
