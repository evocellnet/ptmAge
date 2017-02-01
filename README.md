
# PTM age

Estimating age of PTMs based on a list of sites mapped to Ensembl Compara organisms. For the moment, the pipeline works only in
phosphorylation.


## Dependencies

- Python
- Biopython
- phyml (available in home/linuxbrew)
- newick-utils (available in home/linuxbrew)
- phytools (R package)
- ghostscript
- wget


## Setup

The pipeline takes as input 2 files:

- `data/species_list.txt` a list of species used to generate the phylogenetic trees.
- `data/all_sites_ptmdb.tab` a list of PTMs in some of these species.

Additionally the next variables need to be updated in the `Makefile`.

```makefile

COMPARARELEASE ?= 86
SPECIESTAG ?= species_n10
ANCTHRESHOLD ?= 0.51

```

## Running the pipeline

The pipeline consists in a few steps that are paralellized using LSF's `bsub`. After some of the steps a number of jobs
will be launched and you should not continue to the next step until all of them are finished.

1- Download all Ensembl Compara data, reduce set of organisms and regenerate alignments and trees 

```bash
make prepare
```

2- Runs ancestral reconstructrion and retrieves stats at 2 different windows (0 and +/-3 residues)

```bash
make ancestral
```

This should produce some stats on the different modified alignment columns.

```bash
wc -l results/stats* 
```

3- Collects ancestral states

```bash
make collectStates 
```

4- Compiles all results in summary file

```bash
make collectOrigins 
```

The resulting summary files will be accessible in the `results` directory.

```bash
cat results/all_origins_species_n10_w0_0.51.tab | grep "Homo sapiens" | cut -f9,9 | sort | uniq -c
```


## Authorship ##

This project was developed by Romain Studer. David Ochoa just prepared the Makefile to automate it.

