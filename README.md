[![Build Status](https://travis-ci.org/nickp60/annofilt.svg?branch=master)](https://travis-ci.org/nickp60/annofilt)[![Coverage Status](https://coveralls.io/repos/github/nickp60/annofilt/badge.svg?branch=master)](https://coveralls.io/github/nickp60/annofilt?branch=master)
```
#get a subset of Complete e coli to analyze
Rscript ../riboSeed/scripts/getCompleteGenomeSubset.R /home/nicholas/GitHub/riboSeed/manuscript_results/entropy/assembly_summary.txt ./colis/ Escherichia coli
gunzip ./colis/genomes/*
cd ./colis/genomes/
# for each of them, run prokka to get annotation
counter=0; for i in *.fna; do prokka --outdir ./${i} --prefix coli13 --compliant --genus Escherichia --species coli --cpus 4 ${i}.fna; counter=$(($counter + 1)); done
# I got bored after 11 genomes.  thats enough for a pangenome?  right? right? why are you shaking your head
# build a pangenome
roary -p 4 -f 11complete_colis -e -r -v ./colis/genomes/*/*.gff

# our test organism's annotation
prokka BA000007.2.fasta --outdir sample_prokka --cpus 4

# we also have a mini assembly to test on (see ribSeed toy genome)
prokka --outdir ./assembly_sample --compliant --genus Escherichia --species coli --cpus 4 ../riboSeed/manuscript_results/simulated_genome/test_consensus/final_de_fere_novo_assembly/contigs.fasta
```
![annofilt](https://github.com/nickp60/annofilt/blob/master/icon/icon.svg)

# The Problem
Pangenomes from genome assemblies can be befuddled by missassemblies of genes, expecially those truncated by contig breaks.

# The Solution
`annofilt` is used to filter annotations that appear to be truncated, based on BLAST comparison with a pangenome generated from closed genomes.  Briefly, the algorithm proceeds as follows:


```
for gene in assembly:
   blast against pangenome database
   if the hit passes thresholds set by user:
     retain
   else:
     reject
```


# Building a reference pangenome of trusted genes
To verify the length of annotated genes, we compare annotation length, alignement coverage, and evalue to a pangenme built of well-currated annotations for a given strain.  To build a pangenome for your strain of interest, do the following:
1. Download as many complete genomes (in gff format) from RefSeq as desired (minimum of 10?, maybe?)
2. Run Roary.  This is a good time to explore their stringincy options for percentage identity (which defaults to 95%)
3. Move the `pan_genome_reference.fa` file to a convenient location for use with annofilt.  This contains a representative nucleotide sequences for each gene in the core.

# Running
`annofilt` has three modes:
1 *Normal* (fastest) - check annotations at the beginning and end of contigs
2 *--local_quick* (medium) - use Prokka's protein multifasta  of all genes when blasting; saves the step of writing the genes to disk, but jobs cant be distributed
3. *--full* (slowest) - genes are blasted individually; this gives more control with job hanndling (future versions will hopefully work with SGE, slurm, etc.)

# Quick Start

[Test data can be downloaded here](https://zenodo.org/record/1196324/files/annofilt_test_data_archive.tar.gz)

The test data contains a pangenome of 11 *E. coli* genomes, as well as a complete genome annotated with Prokka, and a toy genome assembly also annotated with Prokka.  To run `annofilt` with the test data, run the following command:

```
annofilt annofilt_test_data_archive/11complete_colis/pan_genome_reference.fa ./annofilt_test_data_archive/assembly_sample/ -o outdir -v 1
```

# Results files
- `all_loci.txt`, `good_loci.txt`, `bad_loci.txt` - these are newline-delimited files containing all locus tags, those passing the user thresholds, and those failing to pass the thresholds, respectively.
- `blast_cmds` - text file contatining BLAST commands used
- `merged_results.tab` - tab-delimitted file containing all blast results before filtering
- `filtered_hits.csv` - comma-delimitted file containing all blast results after filtering
- `*x*.gbk` - GenBank file containing annotations for genes passing thresholds
- `*x*.gff` - gff3 file containing annotations for genes passing thresholds, for use with Roary
