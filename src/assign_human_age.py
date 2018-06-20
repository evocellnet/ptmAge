#!/usr/bin/env python

import os
import sys

gene2treeFile = sys.argv[1]
allsitesFile = sys.argv[2]
species_tag = sys.argv[3]
region = sys.argv[4]
threshold = sys.argv[5].replace(".", "")

# Load family definition
gene_dict = {}
family_list = []
file_in = open(gene2treeFile, "r")
while 1:
    line = file_in.readline()
    if line == "":
        break
    #print line
    line = line.rstrip()
    tab = line.split("\t")
    gene = tab[0]
    path = tab[2]
    gene_dict[gene] = path
    family_list.append(path)
file_in.close()

family_list = list(set(family_list))
family_list.sort()

print "# Number of genes in Compara: ", len(gene_dict)
print "# Number of families in Compara", len(family_list)

#sys.exit()

# Load ubiquinated sites
ubisite_list = []
file_in = open(allsitesFile, "r")
while 1:
    line = file_in.readline()
    if line == "":
        break
    #print line
    tab = line.split("\t")
    gene = tab[1]
    site = tab[2]
    gene_site = gene+"_"+site
    ubisite_list.append(gene)
file_in.close()


ubisite_list= list(set(ubisite_list))
ubisite_list.sort()

print "# Number of genes with (at least 1) ptm sites: ", len(ubisite_list)


# file_out = open(sys.argv[4], "w")
sys.stdout.write("family"+"\t"+"column"+"\t"+"species"+"\t"+ \
                 "gene"+"\t"+"position"+"\t"+"oldest_node"+"\t"+ \
                 "oldest_score"+"\t"+"common_ancestor_tax"+"\t"+ \
                 "common_ancestor_name"+"\n")
for path in family_list:
    # path = gene_dict[gene]
    family = path.split("/")[-1]
    origin_file = path+"/"+species_tag+"/"+"region_w"+region+"/"+family+"_ptm_continue_region_w"+region+"_origins_"+threshold+".txt"
    if os.path.exists(origin_file):
        #print gene, "File exists"
        file_in = open(origin_file, "r")
        while 1:
            line = file_in.readline()
            if line == "":
                break
            line = line.rstrip()
            #print line
            sys.stdout.write(line+"\n")
        file_in.close()
    #else:
    #    print gene, " no file", path+"/"+species_tag+"/"+family+"_ubi_continue_region_w0_origins.txt"
    #    #pass

# file_out.close()
