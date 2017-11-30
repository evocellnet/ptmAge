#!/usr/bin/env python2

import os
import sys
from Bio import SeqIO


# Define species dataset to use
modified_sites_file = sys.argv[1]
ensembl_families_dictfile = sys.argv[2]
species_tag = sys.argv[3]
tree_directory = sys.argv[4]
region = sys.argv[5]
predicted_site = sys.argv[6]

def rgb_to_hex(rgb):
   return '#%02x%02x%02x' % rgb

# Load modified sites per genes
ptm_dict = {}
file_in = open(modified_sites_file, "r")
while 1:
    line = file_in.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t")
    ensembl_id = tab[1]
    site = tab[2]

    # Upload data
    site_list = []
    if ensembl_id in ptm_dict:
        site_list = ptm_dict[ensembl_id]
    site_list.append(site)
    ptm_dict[ensembl_id] = site_list
file_in.close()

# Load modified sites per genes
ptm_predicted_dict = {}
if predicted_site != "":
    file_in = open(predicted_site, "r")
    while 1:
        line = file_in.readline()
        if line == "":
            break
        line = line.rstrip()
        #sys.stderr.write(line + "\n")
        tab = line.split("\t")
        ensembl_id = tab[0]
        site = tab[1]
        proba = tab[4]
    
        # Upload data
        site_dict = {}
        if ensembl_id in ptm_predicted_dict:
            site_dict = ptm_predicted_dict[ensembl_id]
        site_dict[site] = proba
        ptm_predicted_dict[ensembl_id] = site_dict
    file_in.close()

# Load Ensembl families
family_with_ptm = []
ensembl_dict = {}
file_in = open(ensembl_families_dictfile, "r")
while 1:
    line = file_in.readline()
    if line == "":
        break
    tab = line.split()
    ensembl_id = tab[0]
    family = tab[1]

    # Upload data
    ensembl_id_list = []
    if family in ensembl_dict:
        ensembl_id_list = ensembl_dict[family]
    ensembl_id_list.append(ensembl_id)
    ensembl_dict[family] = ensembl_id_list

    # Check if family contains at least 1 ubi sites
    # and add it to the list
    if ensembl_id in ptm_dict:
        family_with_ptm.append(family)
file_in.close()

# Remove duplicated families and sort the list
family_with_ptm = list(set(family_with_ptm))
family_with_ptm.sort()



# Write output
counter = 0

# print "Family"+"\t"+"Column"+"\t"+"Human_genes"+"\t"+"Human_sites"+"\t"+ \
#     "Human_aa"+"\t"+"Human_ubi"+"\t"+"N_genes"+"\t"+"N_aa"+"\t"+"N_Lys"+"\t"+ \
#     "N_ubi_sites"+"\t"+"Ratio ubi_sites / N_aa"+"\t"+"Ratio ubi_sites / N_Lys"+ \
#     "\t"+"All_genes+"\n"


# Write output
# file_out = open(sys.argv[4], "w")
sys.stdout.write("Family"+"\t"+"Column"+"\t"+ \
                 "Human_genes"+"\t"+"Human_site"+"\t"+ \
                 "Human_aa"+"\t"+"Human_phospho"+"\t"+ \
                 "N_genes"+"\t"+"N_aa"+"\t"+"N_STY"+"\t"+ \
                 "N_ubi_sites"+"\t"+"Ratio phosphosites / N_aa"+ \
                 "\t"+"Ratio phosphosites / N_STY"+"\t"+"All_genes"+"\n")

# family_with_ptm = ["ENSTREE_05500"]

for family in family_with_ptm:
    tree_number = int(family.split("_")[1])
    supra_folder = "ENSTREE_"+str(int(str(tree_number).zfill(5)[0:2])+1).zfill(2)+"000"

    # Load sequence file (fasta)
    sequence_dict = {}
    msa_file = tree_directory+"/"+supra_folder+"/"+family+"/"+species_tag+"/"+family+"_reduced.aa.nogap.fasta"
    handle = open(msa_file, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        sequence_dict[record.id] = str(record.seq)
    handle.close()

    column_with_ptm = []
    corresponding_dict = {}
    human_gene = "NA"
    # Screen all sequences
    for ensembl_id in sequence_dict:

        # Define human geen
        if ensembl_id[0:5] == "ENSP0":
            human_gene = ensembl_id

        # Load sequence
        sequence = sequence_dict[ensembl_id]

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

        # Check if ensembl_id contains ubi sites
        if ensembl_id in ptm_dict:
            for site in ptm_dict[ensembl_id]:
                column = corr_dict[int(site)-1]
                if(sequence_dict[ensembl_id][int(column)] in ["S", "T", "Y"]):
                    column_with_ptm.append(str(column)) # Add column to list
                
    # Remove duplicate columns
    column_with_ptm = list(set(column_with_ptm))
    column_with_ptm.sort(key=int)

    scoring_dict = {}
    
    for column_central in column_with_ptm:
        column_central = int(column_central)
        #column = column-1
        n_aa = 0
        n_u_sites = 0

        human_gene_list = []
        human_site_list = []
        human_aa_list = []
        human_ubi_list = []

        aa_list = []
        for ensembl_id in sequence_dict:
            #print column
            # Get main aa
            aa = sequence_dict[ensembl_id][column_central]
            aa_list.append(aa)

            # Define human geen
            human_gene = "NA"
            human_site = "NA"
            human_aa = "NA"
            human_ubi = "NO"
            if ensembl_id[0:5] == "ENSP0":
                human_gene = ensembl_id
                if column in corresponding_dict[ensembl_id]:
                    human_site = str(corresponding_dict[ensembl_id][column]+1)
                human_aa = aa
                if ensembl_id in ptm_dict:
                    if str(human_site) in ptm_dict[ensembl_id]:
                        human_ubi = "YES"
            human_gene_list.append(human_gene)
            human_site_list.append(human_site)
            human_aa_list.append(human_aa)           
            human_ubi_list.append(human_ubi)

                        
            aa_vector = []
            value_vector = []

            # Scan across region
            for column in range(column_central-int(region), column_central+int(region)+1):
                if column in corresponding_dict[ensembl_id]:
                    aa = sequence_dict[ensembl_id][column]
                    aa_vector.append(aa)

                    score = 0.0
                    if aa in ["S", "T", "Y"]:
                        score = 0.5
                    site = str(corresponding_dict[ensembl_id][column]+1)
                    if ensembl_id in ptm_predicted_dict:
                        if site in ptm_predicted_dict[ensembl_id]:
                            score = ptm_predicted_dict[ensembl_id][site]
                    if ensembl_id in ptm_dict:
                        # print ensembl_id, site, ptm_dict[ensembl_id]
                        if str(site) in ptm_dict[ensembl_id]:
                            n_u_sites = n_u_sites+1
                            score = 1.0
                    value_vector.append(score)
                    column_dict = {}

            # Get max score
            best_score = 0
            if len(value_vector) > 0:
                best_score = max(value_vector)
            column_dict = {}
            if ensembl_id in scoring_dict:
                column_dict = scoring_dict[ensembl_id]
            column_dict[int(column_central)] = best_score
            scoring_dict[ensembl_id] = column_dict

            
        n_aa = len(aa_list)
        n_k = aa_list.count("S")+aa_list.count("T")+aa_list.count("Y")
        # print family+"\t"+str(column)
        # print aa_list
        # print family+"\t"+str(column+1)+"\t"+",".join(human_gene_list)+"\t"+",".join(human_site_list)+"\t"+",".join(human_aa_list)
        # print "N aa", n_aa
        # print "N K", n_k
        # print family+"\t"+str(column+1)+"\t"+",".join(human_gene_list)+"\t"+",".join(human_site_list)+"\t"+",".join(human_aa_list)+"\t"+",".join(human_ubi_list)+"\t"+str(len(sequence_dict))+"\t"+str(n_aa)+"\t"+str(n_k)+"\t"+str(n_u_sites)+"\t"+str(round(float(n_u_sites)/n_aa, 2))+"\t"+str(round(float(n_u_sites)/n_k, 2))+"\t"+",".join(sequence_dict.keys())

        sys.stdout.write(family+"\t"+str(column_central+1)+"\t"+";".join(human_gene_list)+"\t"+ \
                         ";".join(human_site_list)+"\t"+";".join(human_aa_list)+"\t"+ \
                         ";".join(human_ubi_list)+"\t"+str(len(sequence_dict))+"\t"+ \
                         str(n_aa)+"\t"+str(n_k)+"\t"+str(n_u_sites)+"\t"+ \
                         str(round(float(n_u_sites)/n_aa, 2))+"\t"+ \
                         str(round(float(n_u_sites)/n_k, 2))+"\t"+ \
                         ";".join(sequence_dict.keys())+"\n")
        #file_out.write(family+"\t"+str(column)+"\t"+human_gene+"\t"+str(human_site)+"\t"+str(human_aa)+"\t"+str(human_ubi)+"\t"+str(len(sequence_dict))+"\t"+str(n_aa)+"\t"+str(n_k)+"\t"+str(n_u_sites)+"\t"+str(round(float(n_u_sites)/n_aa, 2))+"\t"+str(round(float(n_u_sites)/n_k, 2))+"\n")
        #print "ENSTREE_"+ensembl_release+"/"+supra_folder+"/"+family, str(column+1), n_u_sites

    # Write Jalview output
    jalview_out = open(tree_directory+"/"+supra_folder+"/"+family+"/"+species_tag+"/"+family \
                       +".annotations_region_w"+str(region)+".txt", "w")
    jalview_out.write("phospho\tred\n")
    for ensembl_id in sequence_dict:
        for column in column_with_ptm:
            column = int(column)
            #print ensembl_id, column, corresponding_dict[ensembl_id]
            if column in corresponding_dict[ensembl_id]:
                site = corresponding_dict[ensembl_id][column]+1
                score = float(scoring_dict[ensembl_id][column])

                colour = "white"
                gradient = 255
                if score > 0:
                    gradient = int(255*(1-score))
                if gradient > 240:
                    gradient = 240
                colour = rgb_to_hex((gradient,gradient,gradient))
                colour = colour[1:]
                if ensembl_id in ptm_dict:
                    if str(site) in ptm_dict[ensembl_id]:
                        colour = "red"                
                #print ensembl_id, column, score
                if score > 0:
                    jalview_out.write(str(score)+"\t"+ensembl_id+"\t" \
                                      +"-1"+"\t"+ \
                                      str(site)+"\t"+str(site)+"\t"+colour+"\n")
    jalview_out.close()


    #Preparing ancestral reconstruction
    basename = tree_directory+"/"+supra_folder+"/"+family+"/"+species_tag+"/"+"region_w"+str(region)+"/"+family
    os.system("mkdir -p "+tree_directory+"/"+supra_folder+"/"+family+"/"+species_tag+"/"+"region_w"+str(region))
    
    # Write input file for ancestral state reconstruction
    score_file_out = open(basename+"_phospho_continue_region_w"+str(region)+".txt", "w")
    # print "Gene"+"\t"+"\t".join(column_with_ptm)
    real_column_with_ptm = [str(int(column)+1) for column in column_with_ptm]
    score_file_out.write("Gene"+"\t"+"\t".join(real_column_with_ptm)+"\n")
    for ensembl_id in scoring_dict:
        new_line = ensembl_id
        for column in column_with_ptm:
            score = "0.0"
            if int(column) in scoring_dict[ensembl_id]:
                score = scoring_dict[ensembl_id][int(column)]
            new_line = new_line+"\t"+str(score)
        score_file_out.write(new_line+"\n")
    score_file_out.close()

    # Write R script to build ancestral reconstruction
    tree_file = tree_directory+"/"+supra_folder+"/"+family+"/"+species_tag+"/"+family+"_reduced.aa.nogap.ultrametric.tree"
    ubi_file = open(basename+"_phospho_continue_region_w"+str(region)+".r", "w")
    ubi_file.write("if (!require(\"ape\")) {\n")
    ubi_file.write("    library(ape)\n")
    ubi_file.write("}\n")
    ubi_file.write("if (!require(\"phytools\")) {\n")
    ubi_file.write("    library(phytools)\n")
    ubi_file.write("}\n")
    ubi_file.write("phosphotree <- read.tree(\""+tree_file+"\")\n")
    ubi_file.write("phosphotree <- multi2di(phosphotree)\n")
    ubi_file.write("phosphodata <- read.table(\""+family+"_phospho_continue_region_w"+ \
                       str(region)+".txt\",head=T)\n")
    ubi_file.write("phosphodata <- phosphodata[match(phosphotree$tip.label, phosphodata$Gene),]\n")
    for column in column_with_ptm:
        yeast_position = "NA"
        yeast_aa = "NA"
        # if yeast_gene in corr_alignment_reverse:
        #     if int(site) in corr_alignment[yeast_gene]:
        #         yeast_position = corr_alignment[yeast_gene][int(site)]
        #         yeast_aa = gene_dict[yeast_gene][int(site)-1]

        ubi_file.write("char1<-phosphodata$X"+str(int(column)+1)+"\n")
        ubi_file.write("names(char1)<-phosphodata$Gene\n")
        # ubi_file.write("ancestral_reconstruction <- ace(char1, phosphotree, type=\"continuous\", method=\"ML\")\n")
        ubi_file.write("ancestral_reconstruction <- fastAnc(phosphotree, char1, vars=TRUE, CI=TRUE)\n")
        ubi_file.write("phosphotree$node.label <- ancestral_reconstruction$ace\n")
        ubi_file.write("write.tree(phosphotree,\""+family+"_phospho_"+str(int(column)+1).zfill(5)+"_continue_region_w"+str(region)+".tree\")\n")
        ubi_file.write("MLfinal <- c(char1, ancestral_reconstruction$ace)\n")

        # Write PDF
        ubi_file.write("pdf(\""+family+"_phospho_"+str(int(column)+1).zfill(5)+"_continue_region_w"+str(region)+".pdf\", height=1.5+(nrow(phosphodata)/5))\n")
        ubi_file.write("plot(phosphotree, label.offset=0.03, cex=0.6)\n")
        #ubi_file.write("tiplabels(thermo = char1, piecol = c(\"#AACCFF\",\"white\"), cex=0.6, width=.02, height=.5)\n")
        #ubi_file.write("nodelabels(thermo = ancestral_reconstruction$ace, piecol = c(\"#AACCFF\",\"white\"), cex=0.6, width=.02, height=.5)\n")
        ubi_file.write("tiplabels(thermo = char1, piecol = c(\"red\",\"white\"), cex=0.6, width=.02, height=.5)\n")
        ubi_file.write("nodelabels(thermo = ancestral_reconstruction$ace, piecol = c(\"red\",\"white\"), cex=0.6, width=.02, height=.5)\n")
        ubi_file.write("title(main = \"Tree for site "+str(int(column)+1)+" ("+str(yeast_aa)+str(yeast_position)+") in alignment "+family+"\")\n")
        ubi_file.write("dev.off()\n")

    # Close file
    ubi_file.close()

    # Execute script
    os.chdir(tree_directory+"/"+supra_folder+"/"+family+"/"+species_tag+"/"+"region_w"+str(region))
    os.system("bsub -M 4000 -R 'rusage[mem=4000]' -o /dev/null -e /dev/null sh -c \'rm -fr *.pdf;"+ \
              " rm -fr *_continue_region_w"+region+".tree;" + \
              " R CMD BATCH "+family+"_phospho_continue_region_w"+str(region)+".r;"+\
              "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="+ \
              family+"_phospho_continue_region_w"+str(region)+"_trees.pdf "+ \
              family+"_phospho_*_continue_region_w"+str(region)+".pdf\' > /dev/null")
    os.chdir("../../../../../")
    #sys.exit()
# Close file
# file_out.close()
