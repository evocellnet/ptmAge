
# PTM age

Estimating age of PTMs based on a list of sites mapped to Ensembl Compara organisms. For the moment, the pipeline works
only in phosphorylation.


## Dependencies

- Python (v2 & v3)
- Biopython (v2 & v3)
- phyml (available in (home/linux)brew)
- newick-utils (available in (home/linux)brew)
- phytools (R package)
- ghostscript
- wget


## Setup

The pipeline takes as input 2 files:

- `data/species_list.txt` a list of species used to generate the phylogenetic trees.
- `data/all_sites_ptmdb.tab` a list of PTMs in some of these species.

Additionally the next variables need to be updated in the `Makefile`.

```makefile

COMPARARELEASE ?= 86 # Compara release
SPECIESTAG ?= species_n10 #tag to distinguish the organisms list
ANCTHRESHOLD ?= 0.51 # threshold for ancestral reconstruction
SAMPLESIZE ?= 4000 #number of positives and negatives: Total 2x

```

### Running the pipeline ###

The pipeline consists in a few steps that are paralellized using LSF's `bsub`. After some of the steps a number of jobs
will be launched and you should not continue to the next step until all of them are finished.

1- Download all Ensembl data, filter the set of organisms of interest and regenerate alignments and trees 

```bash
make prepare
```

The next few steps calculate the probabilities of all STY to be phosphorylated according to a SVM trained using
netphorest. The default pipeline trains the models using the available data in each organism but it only calculates the
probabilities based on the human model. Some of the steps require some time to run or implemented in the LSF system,
therefore you will have to wait and check that all jobs are finished before runing the next instruction. If you are ok
using 0.5 as a probability for all phosphoacceptor residues you can jump to step X (WARNING: call to `make ancestral`
not yet implemented to run without probabilities).


2- Run Netphorest in every organism

```bash

make netphorest_proteomes
#in LSF "bsub -n 10 lsmake netphorest_proteomes""

```

3- Prepare feature files, training set and test set to train the SVM. The total number of positives and negatives is
contained in the variable `SAMPLESIZE`.

```bash
make prepareAllSVM
```

4- Train and test the SVM based on the previously generated files

```bash
make trainSVM
```

5- Prepare full dataset of STYs to classify

```bash
make prepareFullSVM
```

5- Prepare full dataset of STYs to classify

```bash
make prepareFullSVM
```

6- Classify the full set of sites using the human ST model

```bash
make classifyFull
```

7- Convert the predictions into probabilities

```bash
make getProbabilities
```

If everything worked fine the resulting probabilities based on the SVM can be found in `temp/allSTY_probabilities.tab`.


8- Now we will run the ancestral reconstruction. The next step will retrieve stats at 2 different windows (0 and +/-3
residues).

```bash
make ancestral
```

This step produces some stats on the different modified alignment columns.

```bash
wc -l results/stats* 
```

9- Collects ancestral states

```bash
make collectStates 
```

10- Compiles all results in summary file

```bash
make collectOrigins 
```

The resulting summary files will be accessible in the `results` directory.

```bash
cat results/all_origins_species_n10_w0_0.51.tab | grep "Homo sapiens" | cut -f9,9 | sort | uniq -c
```


## Authorship ##

This project was developed by Romain Studer. David Ochoa just prepared organized the scripts in a common project and
wrote the Makefile to automate the process.

