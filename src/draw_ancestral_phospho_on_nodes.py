#!/usr/bin/env python2

import os
import sys
import re
from Bio import SeqIO
from pytreetools import *

ensemblGeneTree = sys.argv[1]
columnsFile = sys.argv[2]
outputFile = sys.argv[3] #output
anc_threshold = float(sys.argv[4]) #ancestral reconstruction threshold
region = sys.argv[5] #Region window
ncbi_folder = sys.argv[6] # Define taxbrowser folder
allsitesFile = sys.argv[7]

# Define base name
basename = ensemblGeneTree.split(".")[0]


# Definition of the classe Node
class Node:
    """Noeud"""
    def __init__(self):
        self.tax_id = 0       # Number of the tax id.
        self.parent = 0       # Number of the parent of this node
        self.children = []    # List of the children of this node
        self.tip = 0          # Tip=1 if it's a terminal node, 0 if not.
        self.name = ""        # Name of the node: taxa if it's a terminal node, numero if not.       
    def genealogy(self):      # Trace genealogy from root to leaf
        ancestors = []        # Initialise the list of all nodes from root to leaf.
        tax_id = self.tax_id  # Define leaf
        while 1:
            if name_object.has_key(tax_id):
                ancestors.append(tax_id)
                tax_id = name_object[tax_id].parent
            else:
                break
            if tax_id == "1":
                # If it is the root, we reached the end.
                # Add it to the list and break the loop
                ancestors.append(tax_id)
                break
        return ancestors # Return the list

# Function to find common ancestor between two nodes or more
def common_ancestor(taxonomic_tree,node_list):
    name_object = taxonomic_tree
    list1 = name_object[node_list[0]].genealogy()  # Define the whole genealogy of the first node
    for node in node_list:
        list2 = name_object[node].genealogy()      # Define the whole genealogy of the second node
        ancestral_list = []                             
        for i in list1:
            if i in list2:                         # Identify common nodes between the two genealogy
                ancestral_list.append(i)                 
        list1 = ancestral_list                     # Reassing ancestral_list to list 1.
    common_ancestor = ancestral_list[0]            # Finally, the first node of the ancestra_list is the
                                                   # common ancestor of all nodes.
    return common_ancestor                         # Return a node


##################################
#                                #
#   Read NCBI taxonomy files     #
#                                #
##################################

######################
# 
# Load names defintion

name_dict = {}          # Initialise dictionary with TAX_ID:NAME
name_dict_reverse = {}  # Initialise dictionary with NAME:TAX_ID

# Load  NCBI names file ("names.dmp")
name_file =  open(ncbi_folder+"/names.dmp","r")
while 1:
    line = name_file.readline()
    if line == "":
        break
    line = line.rstrip()
    line = line.replace("\t","")
    tab = line.split("|")
    if tab[3] == "scientific name":
        tax_id, name = tab[0], tab[1]     # Assign tax_id and name ...
        name_dict[tax_id] = name          # ... and load them
        name_dict_reverse[name] = tax_id  # ... into dictionaries
    if tab[3] == "synonym":
        tax_id, name = tab[0], tab[1]     # Assign tax_id and name ...
        name_dict_reverse[name] = tax_id  # ... into dictionaries    
name_file.close()

#print name_dict_reverse
#print name_dict_reverse["Homo sapiens"]


# Load all sites
realphospho = {}
file_in = open(allsitesFile, "r")
while 1:
    line = file_in.readline()
    if line == "":
        break
    line = line.rstrip()
    #print line
    tab = line.split("\t")
    gene = tab[1]
    site = tab[2]
    sitelist=[]
    if gene in realphospho:
        sitelist=realphospho[gene]
    sitelist.append(site)
    realphospho[gene]=sitelist    
file_in.close()

    
######################
# 
# Load NCBI taxonomy as tree

# Define taxonomy variable
global name_object
name_object = {}

# Load taxonomy NCBI file ("nodes.dmp")
taxonomy_file = open(ncbi_folder+"/nodes.dmp","r")
while 1:
    line = taxonomy_file.readline()
    if line == "":
        break
    #print line
    line = line.replace("\t","")
    tab = line.split("|")
    
    tax_id = str(tab[0])
    tax_id_parent = str(tab[1])
    division = str(tab[4])

    # Define name of the taxid
    name = "unknown"
    if tax_id in name_dict:
        name = name_dict[tax_id]
    
    if not name_object.has_key(tax_id):
        name_object[tax_id] = Node()
    name_object[tax_id].tax_id   = tax_id        # Assign tax_id
    name_object[tax_id].parent   = tax_id_parent # Assign tax_id parent
    name_object[tax_id].name     = name          # Assign name
    
    if  tax_id_parent in name_object:
        children = name_object[tax_id].children  # If parent is is already in the object
        children.append(tax_id)                  # ...we found its children.
        name_object[tax_id].children = children  # ... so add them to the parent
taxonomy_file.close()

# print "NCBI taxonomy loaded"

# list_all_descendant = []
# list_all_descendant_name = ["Homo sapiens", "Mus musculus"]
# for species in list_all_descendant_name:
#     tax_id = name_dict_reverse[species]
#     list_all_descendant.append(tax_id)
# common_ancestor_all = common_ancestor(name_object, list_all_descendant)
# #print common_ancestor_all
# #print name_dict[common_ancestor_all]

#sys.exit()




# Load alignment

# Load sequence file
sequence_dict = {}
corresponding_dict = {}
msa_file = basename+"_reduced.aa.nogap.fasta"
handle = open(msa_file, "rU")
for record in SeqIO.parse(handle, "fasta"):
    ensembl_id = record.id
    sequence = str(record.seq)
     
    # Define corresponding position
    corr_dict = {}
    corr_dict_rev = {}
    i, j = 0, 0
    for i in range(len(sequence)):
        aa = sequence[i]
        if aa != "-":
            corr_dict[j] = i
            corr_dict_rev[i] = j
            j = j+1
    corresponding_dict[ensembl_id] = corr_dict_rev
handle.close()


# Load taxa name of genes from full Ensembl Gene tree
# gene_to_species_dict = {}
gene_to_taxid = {}
nexus_tree = open(ensemblGeneTree,"r")
while 1:
    line = nexus_tree.readline()
    if line == "" or line[0] == "]":
        break
    # if line[0:3] == "ENS" or "DROSOPHILA" in line or "SACCHAROMYCES" in line or "CAENORHABDITIS" in line:
    if line[0] != "[":
        line = line.rstrip()
        tab = line.split()
        #print tab
        gene = tab[0]
        genus = tab[1].split("_")[0].capitalize()
        species = tab[1].split("_")[1].lower()
        species_name = genus+" "+species
        taxid  = name_dict_reverse[species_name]
        # gene_to_species_dict[gene] = species_name
        gene_to_taxid[gene] = taxid
nexus_tree.close()

# print gene_to_taxid
#print "Ensembl Compara tree loaded"

##############################################
#
# Load column with values in modern

column_list = []
tip_gene_values = {}
file_in = open(columnsFile, "r")
while 1:
    line = file_in.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t")
    if tab[0] == "Gene":
        column_list =  tab[1:]
    else:
        tip_gene = tab[0]
        score_list = tab[1:]
        score_dict = {}
        for i in range(len(column_list)):
            column = column_list[i]
            score = score_list[i]
            #print species, column, score
            score_dict[column] = score
        tip_gene_values[tip_gene] = score_dict
file_in.close()

#print "Scores for modern species loaded"
#print "List of columns to see: ", column_list



##############################################
#
# Iterate all column_list

# print column_list

#sys.exit()

output = open(outputFile, "w")
for column in column_list:
    
    # Load ancestral tree with internal values
    anc_tree_file = "./"+basename+"_phospho_"+str(column).rjust(5, "0")+"_continue_region_w"+region+".tree"
    if os.path.exists(anc_tree_file):
        global ancestral_tree_object
        ancestral_tree_object = {}

        ancestral_tree = ""
        tree_file = open(anc_tree_file,"r")
        while 1:
            line = tree_file.readline() 
            if line == "":
                break
            #print line
            ancestral_tree = line.rstrip()
        tree_file.close()

        ancestral_tree_object = read_tree(ancestral_tree)

        #print "Ancestral tree", ancestral_tree
        # print "Ancestral tree loaded"

        #tab_tree = parsed_tree(ancestral_tree)
        #print tab_tree

        #for node in ancestral_tree_object:
        #    print ancestral_tree_object[node].value()

        tip_list = ancestral_tree_object[1].all_descendant()
        #print "Tip list", tip_list

        for node in tip_list:
            gene = ancestral_tree_object[node].name
            species = name_object[gene_to_taxid[gene]].name
            tip_score = tip_gene_values[gene][column]
            # if tip_score == "1.0":
            if int(column)-1 in corresponding_dict[gene] and gene in realphospho:
                site = corresponding_dict[gene][int(column)-1]+1
                if str(site) in realphospho[gene]:
                    #print "Species, gene, tip_score", species, gene, tip_score
                    #print ancestral_tree_object[node].genealogy()
                    oldest_node = node
                    score = tip_score
                    for internal_node in ancestral_tree_object[node].genealogy():
                        if internal_node in tip_list:
                            #print species_values[species]
                            score = float(tip_gene_values[gene][column])
                        elif ancestral_tree_object[internal_node].bootstrap != "":
                            score =  float(ancestral_tree_object[internal_node].bootstrap)
                        #print internal_node, score
                        if score > anc_threshold and internal_node != 1:
                            oldest_node = internal_node
                            #print "Test", internal_node, score

                    oldest_score = tip_score
                    if ancestral_tree_object[oldest_node].bootstrap != "":
                        oldest_score = round(float(ancestral_tree_object[oldest_node].bootstrap), 2)
                    #node_name = ancestral_tree_object[node]
                    sub_clade_tip = ancestral_tree_object[oldest_node].all_descendant()
                    sub_clade_species = []
                    for tip in sub_clade_tip:
                        sub_clade_species.append(gene_to_taxid[ancestral_tree_object[tip].name])
                    common_ancestor_all = common_ancestor(name_object, sub_clade_species)

                    site = corresponding_dict[gene][int(column)-1]+1
                    #print basename, column, node, oldest_node, ancestral_tree_object[oldest_node].bootstrap, sub_clade_species, common_ancestor_all, name_object[common_ancestor_all].name
                    output.write(basename+"\t"+column+"\t"+species+"\t"+gene+"\t"+str(site)+"\t"+ \
                                 str(oldest_node)+"\t"+str(oldest_score)+"\t"+ \
                                 str(common_ancestor_all)+"\t"+name_object[common_ancestor_all].name+"\n")
output.close()
