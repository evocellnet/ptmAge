#!/usr/bin/env python

import os
import sys

statsfile = sys.argv[1]
species_tag = sys.argv[2]
# ensembl_release = sys.argv[3]
treedir = sys.argv[3]
ncbi_taxonomy = sys.argv[4]
region = sys.argv[5]
ancthreshold = sys.argv[6]
allsitesFile = sys.argv[7]

# Load families
family_list = []
file_in = open(statsfile, "r")
while 1:
    line = file_in.readline()
    if line == "":
        break
    tab = line.split()
    if tab[0] != "Family":
        family = tab[0]
        family_list.append(family)
file_in.close()

family_list = list(set(family_list))
family_list.sort()

# print family_list
# print len(family_list)

for family in family_list:
    tree_number = int(family.split("_")[1])
    supra_folder = "ENSTREE_"+str(int(str(tree_number).zfill(5)[0:2])+1).zfill(2)+"000"
    os.chdir(treedir+"/"+supra_folder+"/"+family+"/"+species_tag+"/region_w"+region)
    os.system("cp ../../"+family+".gnt ./")
    os.system("cp ../"+family+"_reduced.aa.nogap.fasta ./")
    os.system("bsub -e /dev/null -o /dev/null -M 4000 -R 'rusage[mem=4000]' python "+treedir+"/../../src/draw_ancestral_phospho_on_nodes.py "+family+".gnt "+family+"_ptm_continue_region_w"+region+".txt "+family+"_ptm_continue_region_w"+region+"_origins_"+ancthreshold.replace(".","")+".txt "+ancthreshold+" "+region+" "+ncbi_taxonomy+" "+allsitesFile+" > /dev/null")
    # print "bsub -e /dev/null -o /dev/null -M 4000 -R 'rusage[mem=4000]' python "+treedir+"/../../src/draw_ancestral_phospho_on_nodes.py "+family+".gnt "+family+"_ptm_continue_region_w"+region+".txt "+family+"_ptm_continue_region_w"+region+"_origins_"+ancthreshold.replace(".","")+".txt "+ancthreshold+" "+region+" "+ncbi_taxonomy+" "+allsitesFile+" > /dev/null"
    # sys.exit()
