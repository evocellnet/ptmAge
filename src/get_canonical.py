#!/usr/bin/env python3

import os
import sys
import argparse
from Bio import SeqIO

# Parse arguments
parser = argparse.ArgumentParser(prefix_chars='--')
parser.add_argument('--sequences')
parser.add_argument('--compara')
parser.add_argument('--output')

# Extract arguments
args = parser.parse_args()

# 1) Load all proteins from Ensembl proteomes species
compara_sequences_dict = {}
for record in SeqIO.parse(args.compara, "fasta"):
    compara_sequences_dict[record.id] = 1
   
print("Number of sequences loaded from Compara", len(compara_sequences_dict))

# 2) Load all proteins (i.e. from Ensembl proteomes species)
sequences_dict = {}
for record in SeqIO.parse(args.sequences, "fasta"):
    protein_id = record.id
    if record.id[0:3] == "ENS":
        protein_id = ".".join(record.id.split(".")[0:-1])
    sequence = str(record.seq)
    sequences_dict[protein_id] = sequence
   
print("Number of sequences loaded ", len(sequences_dict))
#sys.exit()

# 3) Write sequences
file_out = open(args.output, "w")
for protein_id in sequences_dict:
    if protein_id in compara_sequences_dict:
        file_out.write(">"+protein_id+"\n")
        file_out.write(sequences_dict[protein_id]+"\n")
file_out.close()
