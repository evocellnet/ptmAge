
# Makefile --- 

# Copyright (C) 2016 evolcellnet

# Author: David Ochoa <ochoa@ebi.ac.uk>

# This program is free software, you can redistribute it and/or
# modify it under the terms of the new-style BSD license.

# You should have received a copy of the BSD license along with this
# program. If not, see <http://www.debian.org/misc/bsd.license>.

#################################
# Variables
#################################

COMPARARELEASE ?= 86
SPECIESTAG ?= species_n10
ANCTHRESHOLD ?= 0.51

#Data
SPECIESFILE ?= $(DATADIR)/species_list.txt
ALLSITES ?= $(DATADIR)/all_sites_ptmdb.tab

#################################
# Paths (DIRECTORIES)
#################################

ROOTDIR := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
DATADIR= $(ROOTDIR)/data
TEMPDIR= $(ROOTDIR)/temp
SRCDIR= $(ROOTDIR)/src
RESULTSDIR= $(ROOTDIR)/results

#################################
# Paths (BINS)
#################################

WGET ?= $(shell which wget) -q
PYTHON ?= $(shell which python)
GREP ?= $(shell which grep)
TAR ?= $(shell which tar)
GUNZIP ?= $(shell which gunzip)
RSCRIPT ?= $(shell which Rscript)

COMPARAFTPFASTA ?= ftp://ftp.ensembl.org/pub/release-$(1)/emf/ensembl-compara/homologies/Compara.$(1).protein.aa.fasta.gz
COMPARAFTPTREE ?= ftp://ftp.ensembl.org/pub/release-$(1)/emf/ensembl-compara/homologies/Compara.$(1).protein.nh.emf.gz
NCBITAXONOMYFTP ?= ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

#Temps
COMPARAFASTA ?= $(TEMPDIR)/Compara.$(COMPARARELEASE).protein.aa.fasta
COMPARATREES ?= $(TEMPDIR)/Compara.$(COMPARARELEASE).protein.nh.emf
GENETOTREEFILE ?= $(TEMPDIR)/Gene_to_tree_file.txt
GENETOTREEHUMANFILE ?= $(TEMPDIR)/Gene_to_tree_file_human.txt
GENETREEDIR ?= $(TEMPDIR)/ENSTREE_$(COMPARARELEASE)
NCBITAXONOMY ?= $(TEMPDIR)/taxdump

STATSWINDOW0 ?= $(RESULTSDIR)/stats_$(SPECIESTAG)_w0.tab
STATSWINDOW3 ?= $(RESULTSDIR)/stats_$(SPECIESTAG)_w3.tab

RESULTSWINDOW0 ?= $(RESULTSDIR)/all_origins_$(SPECIESTAG)_w0_0.51.tab
RESULTSWINDOW3 ?= $(RESULTSDIR)/all_origins_$(SPECIESTAG)_w3_0.51.tab

PROTEOMES = $(TEMPDIR)/proteomes
BIOMARTDATASETS = $(TEMPDIR)/biomart_datasets.csv

SIMPLIFYRELEASE = $(shell echo $(1) | sed -e "s/\(.*\)\.p.*/\1/g" | sed -e "s/Equ Cab 2/EquCab2/g" | sed -e "s/JGI 4.2/JGI_4\.2/g")
ENSEMBL_URL = ftp://ftp.ensembl.org/pub/release-$(1)/fasta/$(shell echo $(2) | sed 's/\(.*\) \(.*\)/\L\1_\2/')/pep/$(shell echo $(2) | sed -e "s/\b\(.\)/\u\1/g").$(call SIMPLIFYRELEASE, $(3)).pep.all.fa.gz
ORGSHORTNAME = $(shell echo $(1) | sed 's/\(.\).*[ _]\(.*\)/\L\1\2/')
CSVCUT = $(shell grep $(call ORGSHORTNAME, $(1)) $(BIOMARTDATASETS) | cut -d"," -f$(2) | head -n 1)
ENSEMBL_TARGETS = $(foreach SP,$(SPECIES),ensembl_$(SP))
SPECIES = $(shell sed 's/ /_/g' $(SPECIESFILE))


.PHONY: prepare ancestral collectStates collectOrigins


$(RESULTSDIR):
	mkdir -p $@

$(TEMPDIR):
	mkdir -p $@

#extracting organisms of interest, format conversion, recomputing branch lengths, convert to ultrametric
prepare: | $(TEMPDIR) $(GENETOTREEFILE)
	printf "* Reducing organisms to set of interest\n";\
	$(PYTHON) $(SRCDIR)/ensembl2reduced.py $(GENETOTREEFILE) $(SPECIESFILE) $(SPECIESTAG) 

$(COMPARAFASTA):
	printf "* Downloading Compara Fastas\n";\
	$(WGET) -P $(TEMPDIR) $(call COMPARAFTPFASTA,$(COMPARARELEASE)) -O $@.gz; \
	$(GUNZIP) $@.gz;

$(COMPARATREES):
	printf "* Downloading Compara Tree\n";\
	$(WGET) -P $(TEMPDIR) $(call COMPARAFTPTREE,$(COMPARARELEASE)) -O $@.gz; \
	$(GUNZIP) $@.gz; 

#Extracting alignments and trees from the Ensembl compara files
$(GENETOTREEFILE): $(COMPARATREES) $(COMPARAFASTA)
	printf "* Generating alignments\n";\
	$(PYTHON) $(SRCDIR)/ensembl2aafasta.py $(COMPARAFASTA) $(GENETREEDIR) ;\
	printf "* Generating gene trees\n";\
	$(PYTHON) $(SRCDIR)/ensembl2gnt.py $(COMPARATREES) $(GENETREEDIR) > $@



#Downloads the proteomes for each species from ensembl

ensembl_proteomes:  $(foreach SP,$(SPECIES),ensembl_$(SP))

$(ENSEMBL_TARGETS): ensembl_%: $(PROTEOMES)/%.fasta

#creates proteomes directory
$(PROTEOMES):
	mkdir -p $@

#Downloads the biomart datasets where the genome release is available
$(BIOMARTDATASETS):
	printf "Downloading Biomart datasets...\n";\
	$(RSCRIPT) $(SRCDIR)/biomartDatasets.R ensembl_compara_$(COMPARARELEASE) $@

# Download the Ensembl proteome for each species
$(PROTEOMES)/%.fasta: | $(PROTEOMES)
	printf "Downloading $* from ensembl...\n"
	$(WGET) -P $(PROTEOMES)/$* \
		$(call ENSEMBL_URL,$(COMPARARELEASE),$*,$(call CSVCUT,$*,4)) -O $@.gz
	gunzip $@.gz

#Runs ancestral reconstructrion and retrieves stats at 2 different windows
ancestral: $(STATSWINDOW0) $(STATSWINDOW3)

#Gene to tree mapping in human families
$(GENETOTREEHUMANFILE):
	$(GREP) "^ENSP0" $(GENETOTREEFILE) > $@

#General stats and ancestral reconstruction
$(RESULTSDIR)/stats_$(SPECIESTAG)_w%.tab: $(GENETOTREEHUMANFILE)
	printf "* Extracting stats in conservation\n";\
	$(PHYTON) $(SRCDIR)/compile_data.py $(ALLSITES) $(GENETOTREEHUMANFILE) $(SPECIESTAG) $(GENETREEDIR) $* > $@


#Collect ancestral states
collectStates: states_0 states_3

$(NCBITAXONOMY):
	printf "* Downloading NCBI taxonomy\n";\
	$(WGET) $(NCBITAXONOMYFTP) -O $@.tar.gz; \
	mkdir -p $@ ;\
	$(TAR) -xzf $@.tar.gz -C $@ && rm $@.tar.gz;

states_%: | $(RESULTSDIR) $(NCBITAXONOMY) $(RESULTSDIR)/stats_$(SPECIESTAG)_w%.tab
	$(PYTHON) $(SRCDIR)/collect_ancestral_state.py $(RESULTSDIR)/stats_$(SPECIESTAG)_w$*.tab $(SPECIESTAG) \
	$(GENETREEDIR) $(NCBITAXONOMY) $* $(ANCTHRESHOLD) $(ALLSITES)


#Retrieve origins
collectOrigins: $(RESULTSWINDOW0) $(RESULTSWINDOW3)

$(RESULTSDIR)/all_origins_$(SPECIESTAG)_w%_0.51.tab: | $(RESULTSDIR)
	$(PYTHON) $(SRCDIR)/assign_human_age.py $(GENETOTREEFILE) $(ALLSITES) $(SPECIESTAG) $* $(ANCTHRESHOLD) > $@
