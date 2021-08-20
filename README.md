# Investigating Bacillus subtilis query neighborhoods in metagenome graphs

## Biological questions

### Clouds of species diversity in local environments 

### Prophage insertion sites 



how is this novel? / how compare do existing studies?
- analyze sgc neighborhoods, not assemblies (assemblies will drop too much of the juicy gene bits)
    - if we are working with sgc nbhds, we should probably use protein kmers to get ANI/AAI estimates

    
- MAGs are composite genomes .. how do we handle?
    - must analyze across samples, not w/in a single sample

## Approach

### Clouds of species diversity

+ generate a *B. sub* isolate reference pangenome query
    + Start with a *B. sub* query
    + Run `sourmash prefetch` against GTDBrs202
    + Filter to jaccard similarity of >= 0.1
    + Format for genome-grist
    + genome-grist download prefetch matches
    + charcoal decontaminate genomes 
    + roary generate reference "pangenome" query
+ spacegraphcats query each metagenome CAtlas
+ estimate ANI, AAI, Jaccard
+ MDS plot 
+ pangenome statistics of b. sub
    + multifasta annotation of b. sub catlas (genes present)?
    + protein k-mers?
+ determine additional content recovered by sgc queries
    + megahit assemble query nbhds
    + map nbhd reads back against query nbhds
    + count number of unmapped reads
    + kmer content analysis in full nbhd vs in assemblies
    + ANI estimation from assemblies

## Data

The 605 gut microbiome metagenomes analyzed in this repository were originally analysed in the [2020-ibd](https://github.com/dib-lab/2020-ibd) repository as a meta-cohort of IBD subtypes (CD, UC, and nonIBD).
The download and preprocessing code has been copied over to this workflow, but the pre-processed files have been linked into this repository to save harddisk space.
Specifically, k-mer abundance trimmed pre-processed files have been sym linked into `outputs/abundtrim`, and spacegraphcats metagenome CAtlases have been hard linked into `outputs/sgc_genome_queries`. 
The sample metadata file has been copied into this repository from [here](https://github.com/dib-lab/2020-ibd/blob/master/inputs/working_metadata.tsv).

## Literature

### *B. subtilis*

+ [Genetic Competence Drives Genome Diversity in Bacillus subtilis](https://academic.oup.com/gbe/article/10/1/108/4767717?login=true)
+ [Toward a high-quality pan-genome landscape of Bacillus subtilis by removal of confounding strains](https://academic.oup.com/bib/article-abstract/22/2/1951/5739184)
+ [A pan-genome method to determine core regions of the Bacillus subtilis and Escherichia coli genomes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8156514/)

### Other species:

+ [Genome sequencing of environmental Escherichia coli expands understanding of the ecology and speciation of the model bacterial species](https://www.pnas.org/content/108/17/7200.short)

## Random notes

GTDB extract b sub genomes:

```
grep -i s__Bacillus gtdb-rs202.taxonomy.v2.csv > bacillus.csv
# then add headers, and...
sourmash sig extract --picklist shewanella.csv:ident:ident /group/ctbrowngrp/gtdb/databases/ctb/gtdb-rs202.genomic.k31.zip -o bacillus.zip

```
