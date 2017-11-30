
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
ANCTHRESHOLD ?= 0.46

SAMPLESIZE ?= 4000 #number of positives and negatives: Total 2x

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
PYTHON3 ?= $(shell which python3)
GREP ?= $(shell which grep)
CAT ?= $(shell which cat)
TAR ?= $(shell which tar)
GUNZIP ?= $(shell which gunzip)
RSCRIPT ?= $(shell which Rscript)
NETPHORESTBIN ?= $(SRCDIR)/netphorest
BSUB ?= $(shell which bsub)
SVMLEARN ?= $(SRCDIR)/svm_learn
SVMCLASS ?= $(SRCDIR)/svm_classify

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

RESULTSWINDOW0 ?= $(RESULTSDIR)/all_origins_$(SPECIESTAG)_w0_$(ANCTHRESHOLD).tab
RESULTSWINDOW3 ?= $(RESULTSDIR)/all_origins_$(SPECIESTAG)_w3_$(ANCTHRESHOLD).tab

PROTEOMES = $(TEMPDIR)/proteomes
PROTEOMESCANONICAL = $(TEMPDIR)/proteomes_canonical
NETPHOREST = $(TEMPDIR)/netphorest
BIOMARTDATASETS = $(TEMPDIR)/biomart_datasets.csv
TRAININGSETS = $(TEMPDIR)/training
TESTINGSETS = $(TEMPDIR)/testing
SVMFEATURENAMES = $(TEMPDIR)/feature_names
SVMFEATUREVALUES = $(TEMPDIR)/feature_values
SVMMODELS = $(TEMPDIR)/svmmodels
MODELTESTS = $(TEMPDIR)/svmmodel_tests
TOCLASSIFYSETS = $(TEMPDIR)/toclassify
SVMPREDICTIONS = $(TEMPDIR)/svmpredictions
SVMCLASSSTATS = $(TEMPDIR)/svmclassstats
PROBABILITYMODELS = $(TEMPDIR)/probability_models
SVMPREDICTIONSPROB = $(TEMPDIR)/svmpredictions_probability
ALLPROBABILITIES = $(TEMPDIR)/allSTY_probabilities.tab

SIMPLIFYRELEASE = $(shell echo $(1) | sed -e "s/\(.*\)\.p.*/\1/g" | sed -e "s/Equ Cab 2/EquCab2/g" | sed -e "s/JGI 4.2/JGI_4\.2/g")
ENSEMBL_URL = ftp://ftp.ensembl.org/pub/release-$(1)/fasta/$(shell echo $(2) | sed 's/\(.*\) \(.*\)/\L\1_\2/')/pep/$(shell echo $(2) | sed -e "s/\b\(.\)/\u\1/g").$(call SIMPLIFYRELEASE, $(3)).pep.all.fa.gz
ORGSHORTNAME = $(shell echo $(1) | sed 's/\(.\).*[ _]\(.*\)/\L\1\2/')
CSVCUT = $(shell grep $(call ORGSHORTNAME, $(1)) $(BIOMARTDATASETS) | cut -d"," -f$(2) | head -n 1)
SPECIES = $(shell sed 's/ /_/g' $(SPECIESFILE))
ENSEMBL_TARGETS = $(foreach SP,$(SPECIES),ensembl_$(SP))
ENSEMBL_CANTARGETS = $(foreach SP,$(SPECIES),canonicals_ensembl_$(SP))
NETPHOREST_TARGETS = $(foreach SP,$(SPECIES),netphorest_$(SP))
PREPARESVM_ST_TARGETS = $(foreach SP,$(SPECIES),prepareSVMST_$(SP))
PREPARESVM_Y_TARGETS = $(foreach SP,$(SPECIES),prepareSVMY_$(SP))
TRAIN_ST_PREDICTOR = $(foreach SP,$(SPECIES),train_ST_$(SP))
TRAIN_Y_PREDICTOR = $(foreach SP,$(SPECIES),train_Y_$(SP))
MODELTEST_ST_PREDICTOR = $(foreach SP,$(SPECIES),test_ST_$(SP))
MODELTEST_Y_PREDICTOR = $(foreach SP,$(SPECIES),test_Y_$(SP))
PREPARESVMTOCLASS_ST_TARGETS = $(foreach SP,$(SPECIES),prepareSVMSTtoclass_$(SP))
CLASSIFY_ST_PREDICTOR = $(foreach SP,$(SPECIES),classifyST_$(SP))
PROB_MODELS_TARGETS = $(foreach SP,$(SPECIES),createProbModel_$(SP))
SVM_PROB_TARGETS = $(foreach SP,$(SPECIES),getPredictionProb_$(SP))

.PHONY: prepare ancestral collectStates collectOrigins


$(RESULTSDIR):
	mkdir -p $@

$(TEMPDIR):
	mkdir -p $@

#extracting organisms of interest, format conversion, recomputing branch lengths, convert to ultrametric
#Downloads the proteomes for each species from ensembl
prepare: | $(TEMPDIR) $(GENETOTREEFILE) $(ENSEMBL_TARGETS) $(ENSEMBL_CANTARGETS)
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


$(ENSEMBL_TARGETS): ensembl_%: $(PROTEOMES)/%.fasta
$(ENSEMBL_CANTARGETS): canonicals_ensembl_%: $(PROTEOMESCANONICAL)/%.fasta
$(NETPHOREST_TARGETS): netphorest_%: $(NETPHOREST)/%.netphorest
$(PREPARESVM_ST_TARGETS): prepareSVMST_%: $(TRAININGSETS)/%_ST.train
$(PREPARESVM_Y_TARGETS): prepareSVMY_%: $(TRAININGSETS)/%_Y.train
$(TRAIN_ST_PREDICTOR): train_ST_%: $(SVMMODELS)/%_ST.model 
$(TRAIN_Y_PREDICTOR): train_Y_%: $(SVMMODELS)/%_Y.model
$(MODELTEST_ST_PREDICTOR): test_ST_%: $(MODELTESTS)/%_ST.modeltest 
$(MODELTEST_Y_PREDICTOR): test_Y_%: $(MODELTESTS)/%_Y.modeltest
$(PREPARESVMTOCLASS_ST_TARGETS): prepareSVMSTtoclass_%: $(TOCLASSIFYSETS)/%_ST.full
$(CLASSIFY_ST_PREDICTOR): classifyST_%: $(SVMPREDICTIONS)/%_ST.predictions
$(PROB_MODELS_TARGETS): createProbModel_%: $(PROBABILITYMODELS)/%_ST.probabilitymodel
$(SVM_PROB_TARGETS): getPredictionProb_%: $(SVMPREDICTIONSPROB)/%_ST.predictions

ensembl_proteomes:  $(foreach SP,$(SPECIES),ensembl_$(SP))
ensembl_canonical_proteomes:  $(foreach SP,$(SPECIES),canonicals_ensembl_$(SP))

#they can run with lsmake
netphorest_proteomes: $(foreach SP,$(SPECIES),netphorest_$(SP))

#All this jobs require bsub
prepareAllSVM: $(foreach SP,$(SPECIES),prepareSVMST_$(SP)) $(foreach SP,$(SPECIES),prepareSVMY_$(SP))

trainSVM: $(foreach SP,$(SPECIES),test_ST_$(SP)) $(foreach SP,$(SPECIES),test_Y_$(SP))

#All this jobs require bsub
prepareFullSVM: $(foreach SP,$(SPECIES),prepareSVMSTtoclass_$(SP))

#All this jobs require bsub
classifyFull: $(foreach SP,$(SPECIES),classifyST_$(SP))

#Get probabilities
getProbabilities: $(ALLPROBABILITIES)

#Preparing the training sets
$(TRAININGSETS)/%_ST.train: | $(TRAININGSETS) $(TESTINGSETS) $(SVMFEATURENAMES) $(SVMFEATUREVALUES) $(NETPHOREST)/%.netphorest
	$(BSUB) -o /dev/null -M 25000 -R "rusage[mem=25000]" \
	$(PYTHON3) $(SRCDIR)/netphorest_to_svm.py \
	"--sequences="$(PROTEOMESCANONICAL)/$*.fasta \
	"--psites="$(ALLSITES) \
	"--netphorest="$(NETPHOREST)/$*.netphorest \
	"--species="$* \
	"--out_train="$@ \
	"--out_test="$(TESTINGSETS)/$*_ST.test \
	"--features_names="$(SVMFEATURENAMES)/$*_ST \
	"--features_values="$(SVMFEATUREVALUES)/$*_ST \
	"--sample_size="$(SAMPLESIZE) \
	"--aa_target="ST

$(TRAININGSETS)/%_Y.train: | $(TRAININGSETS) $(TESTINGSETS) $(SVMFEATURENAMES) $(SVMFEATUREVALUES) $(NETPHOREST)/%.netphorest
	$(BSUB) -o /dev/null -M 25000 -R "rusage[mem=25000]" \
	$(PYTHON3) $(SRCDIR)/netphorest_to_svm.py \
	"--sequences="$(PROTEOMESCANONICAL)/$*.fasta \
	"--psites="$(ALLSITES) \
	"--netphorest="$(NETPHOREST)/$*.netphorest \
	"--species="$* \
	"--out_train="$@ \
	"--out_test="$(TESTINGSETS)/$*_Y.test \
	"--features_names="$(SVMFEATURENAMES)/$*_Y \
	"--features_values="$(SVMFEATUREVALUES)/$*_Y \
	"--sample_size="$(SAMPLESIZE) \
	"--aa_target="Y

#Learning models
$(SVMMODELS)/%_ST.model: | $(SVMMODELS) $(TRAININGSETS)/%_ST.train
	$(SVMLEARN) $(TRAININGSETS)/$*_ST.train $@

$(SVMMODELS)/%_Y.model: | $(SVMMODELS) $(TRAININGSETS)/%_Y.train
	$(SVMLEARN) $(TRAININGSETS)/$*_Y.train $@

#Test the model
$(MODELTESTS)/%_ST.modeltest: | $(MODELTESTS) $(SVMMODELS)/%_ST.model
	- $(SVMCLASS) \
	$(TESTINGSETS)/$*_ST.test \
	$(SVMMODELS)/$*_ST.model \
	$@

$(MODELTESTS)/%_Y.modeltest: | $(MODELTESTS) $(SVMMODELS)/%_Y.model
	- $(SVMCLASS) \
	$(TESTINGSETS)/$*_Y.test \
	$(SVMMODELS)/$*_Y.model \
	$@

#Preparing the full sets to classify (only STs)
$(TOCLASSIFYSETS)/%_ST.full: | $(TOCLASSIFYSETS) $(TESTINGSETS) $(SVMFEATURENAMES) $(SVMFEATUREVALUES) $(NETPHOREST)/%.netphorest
	$(BSUB) -o /dev/null -M 20000 -R "rusage[mem=20000]" \
	$(PYTHON3) $(SRCDIR)/netphorest_to_svm.py \
	"--sequences="$(PROTEOMESCANONICAL)/$*.fasta \
	"--psites="$(ALLSITES) \
	"--netphorest="$(NETPHOREST)/$*.netphorest \
	"--species="$* \
	"--out_train=/dev/null" \
	"--out_test="$@ \
	"--features_names=/dev/null" \
	"--features_values=/dev/null" \
	"--sample_size=-1" \
	"--aa_target="ST

$(SVMPREDICTIONS)/%_ST.predictions: | $(SVMPREDICTIONS) $(SVMCLASSSTATS)
	$(BSUB) -o /dev/null \
	$(SVMCLASS) \
	$(TOCLASSIFYSETS)/$*_ST.full \
	$(SVMMODELS)/homo_sapiens_ST.model \
	$@ \
	> $(SVMCLASSSTATS)/$*_ST.stats

$(PROBABILITYMODELS)/%_ST.probabilitymodel: | $(PROBABILITYMODELS)
	$(PYTHON3) $(SRCDIR)/convert_svm_to_proba_extract_AB.py \
	$(MODELTESTS)/$*_ST.modeltest \
	$(TESTINGSETS)/$*_ST.test \
	$@

$(SVMPREDICTIONSPROB)/%_ST.predictions: | $(SVMPREDICTIONSPROB) $(PROBABILITYMODELS)/homo_sapiens_ST.probabilitymodel
	$(PYTHON3) $(SRCDIR)/convert_svm_to_proba_apply_AB.py \
	$(SVMPREDICTIONS)/$*_ST.predictions \
	$(TOCLASSIFYSETS)/$*_ST.full \
	$@ \
	$(PROBABILITYMODELS)/homo_sapiens_ST.probabilitymodel

$(ALLPROBABILITIES): $(foreach SP,$(SPECIES),getPredictionProb_$(SP))
	$(CAT) $(SVMPREDICTIONSPROB)/* > $@

#creates proteomes directory
$(PROTEOMES):
	mkdir -p $@

$(PROTEOMESCANONICAL):
	mkdir -p $@

$(NETPHOREST):
	mkdir -p $@

$(TRAININGSETS):
	mkdir -p $@

$(TESTINGSETS):
	mkdir -p $@

$(SVMFEATURENAMES):
	mkdir -p $@

$(SVMFEATUREVALUES):
	mkdir -p $@

$(SVMMODELS):
	mkdir -p $@

$(TOCLASSIFYSETS):
	mkdir -p $@

$(SVMPREDICTIONS):
	mkdir -p $@

$(SVMCLASSSTATS):
	mkdir -p $@

$(PROBABILITYMODELS):
	mkdir -p $@

$(SVMPREDICTIONSPROB):
	mkdir -p $@

$(MODELTESTS):
	mkdir -p $@

#Downloads the biomart datasets where the genome release is available
# Requires biomaRt and RMySQL
$(BIOMARTDATASETS):
	printf "* Downloading Biomart datasets...\n" ;\
	$(RSCRIPT) $(SRCDIR)/biomartDatasets.R ensembl_compara_$(COMPARARELEASE) $@

# Download the Ensembl proteome for each species
$(PROTEOMES)/%.fasta: | $(PROTEOMES) $(BIOMARTDATASETS)
	printf "* Downloading $* from ensembl...\n" ;\
	$(WGET) -P $(PROTEOMES)/$* \
		$(call ENSEMBL_URL,$(COMPARARELEASE),$*,$(call CSVCUT,$*,4)) -O $@.gz
	gunzip $@.gz

# Requires biopython
$(PROTEOMESCANONICAL)/%.fasta: | $(PROTEOMESCANONICAL) $(PROTEOMES)/%.fasta
	printf "Filtering canonicals for $* from ensembl...\n" ;\
	$(PYTHON3) $(SRCDIR)/get_canonical.py \
	"--sequences="$(PROTEOMES)/$*.fasta \
	"--compara="$(COMPARAFASTA) \
	"--output="$@

$(NETPHOREST)/%.netphorest: | $(NETPHOREST) $(PROTEOMESCANONICAL)/%.fasta
	printf "* Netphorest for $* from ensembl...\n" ;\
	cat $(PROTEOMESCANONICAL)/$*.fasta | $(NETPHORESTBIN) > $@


#Runs ancestral reconstructrion and retrieves stats at 2 different windows
ancestral: $(STATSWINDOW0) $(STATSWINDOW3)

#Gene to tree mapping in human families
$(GENETOTREEHUMANFILE):
	$(GREP) "^ENSP0" $(GENETOTREEFILE) > $@

#General stats and ancestral reconstruction
$(RESULTSDIR)/stats_$(SPECIESTAG)_w%.tab: $(GENETOTREEHUMANFILE)
	printf "* Extracting stats in conservation\n";\
	$(PHYTON) $(SRCDIR)/compile_data.py $(ALLSITES) $(GENETOTREEHUMANFILE) $(SPECIESTAG) $(GENETREEDIR) $* \
	$(ALLPROBABILITIES) > $@

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

$(RESULTSDIR)/all_origins_$(SPECIESTAG)_w%_$(ANCTHRESHOLD).tab: | $(RESULTSDIR)
	$(PYTHON) $(SRCDIR)/assign_human_age.py $(GENETOTREEFILE) $(ALLSITES) $(SPECIESTAG) $* $(ANCTHRESHOLD) > $@


# cat /nfs/research2/beltrao/ochoa/ptmAge/temp/svmpredictions_probability/* > /nfs/research2/beltrao/ochoa/ptmAge/temp/allSTY_probabilities.tab
