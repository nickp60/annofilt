[![Build Status](https://travis-ci.org/nickp60/annofilt.svg?branch=master)](https://travis-ci.org/nickp60/annofilt)[![Coverage Status](https://coveralls.io/repos/github/nickp60/annofilt/badge.svg?branch=master)](https://coveralls.io/github/nickp60/annofilt?branch=master)
<!-- ``` -->
<!-- #get a subset of Complete e coli to analyze -->
<!-- Rscript ../riboSeed/scripts/getCompleteGenomeSubset.R /home/nicholas/GitHub/riboSeed/manuscript_results/entropy/assembly_summary.txt ./colis/ Escherichia coli -->
<!-- gunzip ./colis/genomes/* -->
<!-- cd ./colis/genomes/ -->
<!-- # for each of them, run prokka to get annotation -->
<!-- counter=0; for i in *.fna; do prokka --outdir ./${i} --prefix coli13 --compliant --genus Escherichia --species coli --cpus 4 ${i}.fna; counter=$(($counter + 1)); done -->
<!-- # I got bored after 11 genomes.  thats enough for a pangenome?  right? right? why are you shaking your head -->
<!-- # build a pangenome -->
<!-- roary -p 4 -f 11complete_colis -e -r -v ./colis/genomes/*/*.gff -->

<!-- # our test organism's annotation -->
<!-- prokka BA000007.2.fasta --outdir sample_prokka --cpus 4 -->

<!-- # we also have a mini assembly to test on (see riboSeed toy genome) -->
<!-- prokka --outdir ./assembly_sample --compliant --genus Escherichia --species coli --cpus 4 ../riboSeed/manuscript_results/simulated_genome/test_consensus/final_de_fere_novo_assembly/contigs.fasta -->
<!-- ``` -->
![annofilt](https://github.com/nickp60/annofilt/blob/master/docs/icon/icon.svg)

# The Problem
Pangenomes from genome assemblies can be befuddled by missassemblies of genes. Often the repeats from  multicopy genes cause regions ttat are impossible to assemble from short reads; this results in truncated genes often being found on contig ends.

# The Solution
`annofilt` is used to filter annotations that appear to be truncated, based on BLAST comparison with a pangenome generated from closed genomes.  Briefly, the algorithm proceeds as follows:


```
for gene in assembly:
   blast against pangenome database
   if no hits
     if stringent:
	   reject
     else:
	   retain
   else if hit passes thresholds set by user:
     retain
   else:
     reject
create filtered .gbk file
create filtered .gff file
```


# Building a reference pangenome of trusted genes
To verify the length of annotated genes, we compare annotation length, alignement coverage, and evalue to a pangenme built of well-currated annotations for a given strain.  To build a pangenome for your strain of interest, do the following:

1. Download as many complete genomes from RefSeq as desired (minimum of 10?, maybe?) with `get_compete_genomes``
2. Create pangenome with `make_annofilt_pangenome`.  This is a good time to explore their stringincy options for percentage identity (which defaults to 95%).  If you want, you can adjust the default params using the `--add_roary_args` command.
4. Take a look at the resulting `summary_statistics.txt` file, to make sure nothing looks amiss.
3. Move the `pan_genome_reference.fa` file to a convenient location for use with annofilt.  This contains a representative nucleotide sequences for each gene in the core.

# Installation
```
conda create -n annofilt -c conda-forge -c bioconda prokka roary blast
conda activate annofilt
pip install annofilt
```

# Quick Start

Download 10 random complete genomes
```
get_complete_genomes --genus Escherichia --species coli -o coli_genomes -n 25
```

Annotate them all with Prokka, and then generate a pangenome with Roary
```
make_annofilt_pangenome --genomes coli_genomes/  --output pangenome --threads 4
```

Filter the annotations of an assembly based on this pangenome of trusted genes:

```
annofilt pangenome/pan_genome_reference.fa ./path/to/some/contigs.fasta -o annofilt_results
```

# Running
`annofilt` has three modes:
1. *Normal* (fastest) - check annotations at the beginning and end of contigs
2. *--local_quick* (medium) - use Prokka's protein multifasta  of all genes when blasting; saves the step of writing the genes to disk, but jobs cant be distributed
3. *--full* (slowest) - genes are blasted individually; this gives more control with job hanndling (future versions will hopefully work with SGE, slurm, etc.)

# Quick Start, with sample data

[Test data can be downloaded here](https://zenodo.org/record/1196324/files/annofilt_test_data_archive.tar.gz)

The test data contains a pangenome of 11 *E. coli* genomes, as well as a complete genome annotated with Prokka, and a toy genome assembly also annotated with Prokka.  To run `annofilt` with the test data, run the following command:

```
annofilt annofilt_test_data_archive/11complete_colis/pan_genome_reference.fa ./annofilt_test_data_archive/assembly_sample/ -o outdir -v 1
```

# Results files
- `all_loci.txt`, `good_loci.txt`, `bad_loci.txt` - these are newline-delimited files containing all locus tags, those passing the user thresholds, and those failing to pass the thresholds, respectively
- `nohit_loci.txt` - this newline-delimited file contains genes that failed to get any blast hits, indicating they are not in the core genome
- `blast_cmds` - text file contatining BLAST commands used
- `merged_results.tab` - tab-delimitted file containing all blast results before filtering
- `filtered_hits.csv` - comma-delimitted file containing all blast results after filtering
- `*x*.gbk` - GenBank file containing annotations for genes passing thresholds
- `*x*.gff` - gff3 file containing annotations for genes passing thresholds, for use with Roary


# So what does it do to my assemblies?
I used a subset of the Enterobase E coli collection, where I downloaded a representative from each Acktman sequence types (~1100 strains).

By default, annofilt checks the annotations at the end of each contig. The figure below shows the number of genes searched (2 * number of contigs) in gray, and the number of genes retained is in red.

![searchedvkept](https://github.com/nickp60/annofilt/blob/master/docs/readme_figs/ent2.png)


Here we show the percentage of the searched genes that annofilt retains:

![percsearchedkept](https://github.com/nickp60/annofilt/blob/master/docs/readme_figs/ents.png)

At a more granular scale, here are two instances in mauive alignmens where genes were removed because their identities were ambiguous.  We aligned 3 reference genomes to both the filtered and the unfiltered assemblies for the isolate we're interested in.

![weird_gene3](https://github.com/nickp60/annofilt/blob/master/docs/readme_figs/weird_gene3.png)
See how the genome on top has a large gene where the gene at the end of the contig on bottom is truncated?  The 4th genome (the annofilt one) shows that this partial gene annotation has been removed.

Heres another one:
![weird_gene2](https://github.com/nickp60/annofilt/blob/master/docs/readme_figs/weird_gene2.png)
There is poor homology in that neon green area; the top genome has a syntenous region, but the gene at the end of the contig appears to be truncated.

Overall, in the pangenome we generated with and without annofilt, we reduced the cloud genes by ~5000, and increased the core genome by 70 genes.

![summary_fig](https://github.com/nickp60/annofilt/blob/master/docs/readme_figs/summary.png)


(This comparison also included ~150 of strains of interest to our group)


## Notes for running with Docker
To keep the image size rom being outrageously large, we did not include Prokka in the image.  I maintain a separate Prokka image, which can be obtained from docker hub.  So, we really only recoomend using Docker to run the main annofilt procedure, not using it to run `get_complete_genomes` or `make_annofilt_pangenome`. .
