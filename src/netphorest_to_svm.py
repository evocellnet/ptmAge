#!/usr/bin/env python3

import sys
import argparse
import gzip
import random
from Bio import SeqIO

not_found = 0

# Parse arguments
parser = argparse.ArgumentParser(prefix_chars='--')
parser.add_argument('--sequences')
parser.add_argument('--psites')
parser.add_argument('--netphorest')
parser.add_argument('--out_train')
parser.add_argument('--out_test')
parser.add_argument('--sample_size')
parser.add_argument('--aa_target')
parser.add_argument('--features_values')
parser.add_argument('--features_names')
parser.add_argument('--species')

# Extract arguments
args = parser.parse_args()

# 1) Load all proteins (i.e. from Ensembl proteomes species)
sequences_dict = {}
for record in SeqIO.parse(args.sequences, "fasta"):
    protein_id = ".".join(record.id.split(".")[0:-1])
    sequence = str(record.seq)
    sequences_dict[record.id] = sequence

print("Number of sequences loaded ", len(sequences_dict))
#sys.exit()


# Load phosphosites
site_problem = []
phosphosites_dict = {}
phosphosites_file = open(args.psites, "r")
while 1:
    line = phosphosites_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t")
    #print(tab)

    if tab[0] != "protein_id" and len(tab) > 2:
        species = tab[0]
        gene = tab[1]
        site = tab[2]
        #print(species.replace(" ", "_").lower())
        if species.replace(" ", "_").lower() == args.species:
            if gene in sequences_dict:
                sequence = sequences_dict[gene]
                aa = sequence[int(site)-1]
                if args.aa_target == "ST":
                    if aa in ["S", "T"]:
                        list_of_site = {}
                        if gene in phosphosites_dict:
                            list_of_site = phosphosites_dict[gene]
                        list_of_site[site] = aa
                        #print gene, site, aa
                        phosphosites_dict[gene] = list_of_site
                elif args.aa_target == "Y":
                    if aa in ["Y"]:
                        list_of_site = {}
                        if gene in phosphosites_dict:
                            list_of_site = phosphosites_dict[gene]
                        list_of_site[site] = aa
                        #print gene, site, aa
                        phosphosites_dict[gene] = list_of_site
phosphosites_file.close()


i = 0
for item in phosphosites_dict.items():
    #print item, len(item[1])
    i = i+len(item[1])

#print(phosphosites_dict)
#print(not_found)
print("Number of phosphosites (only "+args.aa_target+"):", i)
print("Number of genes with phosphosites:", len(phosphosites_dict))

#sys.exit()


phosphosites_others_dict = {}
#sys.exit()



# Load phosphosite predictions from Netphorest
kinase_tested_list = []
phosphosite_exp = {}
phosphopred_dict = {}
phosphopred_file = open(args.netphorest, "r")
iterator = 0
while 1:
    line = phosphopred_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split()
    #print(line)
    #print(tab)
    iterator = iterator+1
    # if iterator > 100000:
    #     break
    if len(tab) > 1:
        gene = tab[0]
        #print(str(gene))
        # Check if gene is in phosphosites_dict, to ensure
        # that the pos and neg are from the same dataset
        if (gene in phosphosites_dict) or int(args.sample_size) == -1:
            site = tab[1]
            aa = tab[2]
            gene_site = gene+"_"+site
            predictor = tab[4]
            #if predictor in ["scansite_TurkLabOPL", "nn"]:
            if 1:
                kinase = "_".join(tab[4:8])
                #print(kinase)
                posterior = tab[8]
            if args.aa_target == "ST":
                if aa in ["S", "T"]:
                    if gene in phosphosites_dict:
                        if site in phosphosites_dict[gene]:
                            phosphosite_exp[gene_site] = "+1"
                        else:
                            phosphosite_exp[gene_site] = "-1"
                    elif gene in phosphosites_others_dict:
                        if site not in phosphosites_others_dict[gene]:
                            phosphosite_exp[gene_site] = "-1"
                        else:
                            phosphosite_exp[gene_site] = "-1"
                    else:
                        phosphosite_exp[gene_site] = "-1"

                    if gene_site in phosphosite_exp:
                        kinase_dict = {}
                        if gene_site in phosphopred_dict:
                            kinase_dict = phosphopred_dict[gene_site]
                        kinase_dict[kinase] = posterior
                        phosphopred_dict[gene_site] = kinase_dict

                        # Add kinase to the global list
                    kinase_tested_list.append(kinase)
            if args.aa_target == "Y":
                if aa in ["Y"]:
                    if site in phosphosites_dict[gene]:
                        phosphosite_exp[gene_site] = "+1"
                    elif gene in phosphosites_others_dict:
                        if site not in phosphosites_others_dict[gene]:
                            phosphosite_exp[gene_site] = "-1"
                    else:
                        phosphosite_exp[gene_site] = "-1"


                    if gene_site in phosphosite_exp:
                        kinase_dict = {}
                        if gene_site in phosphopred_dict:
                            kinase_dict = phosphopred_dict[gene_site]
                        kinase_dict[kinase] = posterior
                        phosphopred_dict[gene_site] = kinase_dict

                        # Add kinase to the global list
                    kinase_tested_list.append(kinase)


phosphopred_file.close()

print("Number of prediction from Netphorest:", len(phosphosite_exp))

#print kinase_tested_list
print("N kinase groups ", len(kinase_tested_list))

kinase_tested_list = list(set(kinase_tested_list))
kinase_tested_list.sort()
print("N kinase groups unique ", len(kinase_tested_list))

#sys.exit()


print(len(kinase_tested_list))
kinase_tested_list = list(set(kinase_tested_list))
kinase_tested_list.sort()
#print kinase_tested_list
print(len(kinase_tested_list))

#sys.exit()
all_values = open(args.features_values, "w")
positive_set = []
negative_set = []

count_i = 0
for gene_site in phosphopred_dict:
    #print gene_site

    count_i = count_i+1
    # print(count_i)

    sign = phosphosite_exp[gene_site]
    line = sign
    for i, kinase in enumerate(kinase_tested_list):
        #kinase = kinase_tested_list[i]
        kinase_index = i+1
        kinase_score = 0
        if kinase in phosphopred_dict[gene_site]:
            kinase_score = phosphopred_dict[gene_site][kinase]
        line = line+" "+str(kinase_index)+":"+str(kinase_score)
    line = line+" #"+gene_site
    all_values.write(line+"\n")
    #print line
    if sign == "+1":
        positive_set.append(line)
    elif sign == "-1":
        negative_set.append(line)
all_values.close()

print(len(kinase_tested_list))
print("Total number of positive", len(positive_set))
print("Total number of negative", len(negative_set))

# Write vector list
vector_file = open(args.features_names, "w")
for i, kinase in enumerate(kinase_tested_list):
    kinase = kinase_tested_list[i]
    kinase_index = i+1
    vector_file.write(str(kinase_index)+"\t"+kinase+"\n")
vector_file.close()




sampling_size = int(args.sample_size)

if sampling_size == -1:
    training_sampling_pos = len(positive_set)
    training_sampling_neg = len(negative_set)
    testing_sampling_pos = len(positive_set)
    testing_sampling_neg = len(negative_set)
else:
    training_sampling_pos = int(sampling_size)
    training_sampling_neg = int(sampling_size)
    testing_sampling_pos = int(sampling_size/2)
    testing_sampling_neg = int(sampling_size/2)



# Write svm_output training set
#for ratio in [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,5,7,10]:

svm_train_set = open(args.out_train, "w")
positive_subset = random.sample(positive_set, training_sampling_pos)
for line in positive_subset:
    svm_train_set.write(line+"\n")
negative_subset = random.sample(negative_set, training_sampling_neg)
for line in negative_subset:
    svm_train_set.write(line+"\n")
svm_train_set.close()

print("Number of positive sites used:", len(positive_subset))
print("Number of negative sites used:", len(negative_subset))


print("positive_set", len(positive_set))
print("testing_sampling_pos", testing_sampling_pos)

# Write svm_output test set
svm_test_set = open(args.out_test, "w")
if sampling_size != -1:
    positive_set = list(set(positive_set)-set(positive_subset))
positive_subset = random.sample(positive_set, testing_sampling_pos)
for line in positive_subset:
    svm_test_set.write(line+"\n")
if sampling_size != -1:
    negative_set = list(set(negative_set)-set(negative_subset))
negative_subset = random.sample(negative_set, testing_sampling_neg)
for line in negative_subset:
    svm_test_set.write(line+"\n")
svm_test_set.close()
