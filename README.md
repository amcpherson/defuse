[TOC]

# deFuse

deFuse is a software package for gene fusion discovery using RNA-Seq data. The software uses clusters of discordant paired end alignments to inform a split read alignment analysis for finding fusion boundaries. The software also employs a number of heuristic filters in an attempt to reduce the number of false positives and produces a fully annotated output for each predicted fusion. 

## Publications

The deFuse algorithm and results from an application to ovarian tumours and sarcomas have been published in PLoS Computational Biology: 

[deFuse: An Algorithm for Gene Fusion Discovery in Tumor RNA-Seq Data](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1001138)

deFuse has been used to discover gene fusions in tumour samples for the following publications: 

[Genetic alterations activating kinase and cytokine receptor signaling in high-risk acute lymphoblastic leukemia](http://www.ncbi.nlm.nih.gov/pubmed/22897847)

[14-3-3 fusion oncogenes in high-grade endometrial stromal sarcoma](http://www.ncbi.nlm.nih.gov/pubmed/22223660)

[Whole-genome sequencing identifies genetic alterations in pediatric low-grade gliomas](http://www.ncbi.nlm.nih.gov/pubmed/23583981)

[Recurrent somatic alterations of FGFR1 and NTRK2 in pilocytic astrocytoma](http://www.ncbi.nlm.nih.gov/pubmed/23817572)

[Integrated genome and transcriptome sequencing identifies a novel form of hybrid and aggressive prostate cancer](http://www.ncbi.nlm.nih.gov/pubmed/22294438)

[Poly-gene fusion transcripts and chromothripsis in prostate cancer](http://www.ncbi.nlm.nih.gov/pubmed/22927308)

[Identification of recurrent FGFR3 fusion genes in lung cancer through kinome-centred RNA sequencing](http://www.ncbi.nlm.nih.gov/pubmed/23661334)

[CUX1 is a haploinsufficient tumor suppressor gene on chromosome 7 frequently inactivated in acute myeloid leukemia](http://www.ncbi.nlm.nih.gov/pubmed/23212519)

[Characterization of novel genomic alterations and therapeutic approaches using acute megakaryoblastic leukemia xenograft models](http://www.ncbi.nlm.nih.gov/pubmed/23045605)

[TBL1XR1/TP63: a novel recurrent gene fusion in B-cell non-Hodgkin lymphoma](http://www.ncbi.nlm.nih.gov/pubmed/22496164)

[CBFA2T3-GLIS2 fusion transcript is a novel common feature in pediatric, cytogenetically normal AML, not restricted to FAB M7 subtype](http://www.ncbi.nlm.nih.gov/pubmed/23407549)

[Recurrent mutations in epigenetic regulators, RHOA and FYN kinase in peripheral T cell lymphomas](http://www.ncbi.nlm.nih.gov/pubmed/24413734)

[MHC class II transactivator CIITA is a recurrent gene fusion partner in lymphoid cancers, Nature 2011](http://www.nature.com/nature/journal/v471/n7338/full/nature09754.html)

## Issues

If you are having trouble running deFuse, please open an issue on the [bitbucket issue tracker](https://bitbucket.org/dranew/defuse/issues).

## RNA-Seq Datasets

Four of the granulosa cell tumour datasets analysed in the PLoS Comp. Bio. paper have been submitted to the European Genome-Phenome Archive and can be accessed [here](http://www.ebi.ac.uk/ega/datasets/EGAD00000000046). A mapping between IDs used in the PLoS Comp. Bio. paper and the [New England Journal of Medicine paper](http://www.nejm.org/doi/full/10.1056/NEJMoa0902542) is provided below. 

| N. Engl. J. Med ID | PLoS Comp. Bio. ID |
| ------------------ | ------------------ |
| GCT0026            | GRC4               |
| GCT0028            | GRC1               |
| GCT0077            | GRC2               |
| GCT0078            | GRC3               |

## Setup

To run the deFuse pipeline, you will be required to install some prerequisite software, modify a configuration file, and run a setup script.

### deFuse code

Download and untar/gunzip the deFuse code package to a directory. 

To build the deFuse toolset you must have the **boost c++ development libraries** installed. If they are not installed on your system you can download them from the [boost website](http://www.boost.org/). A full install of boost is not required. The easiest thing to do is to download the latest boost source tar.gz, extract it, then add the extracted path to the CPLUS_INCLUDE_PATH environment variable (in bash, `export CPLUS_INCLUDE_PATH=/boost/directory/:$CPLUS_INCLUDE_PATH`) 

To build the deFuse toolset, cd to the tools subdirectory and type `make`. 

Create a directory for storing the reference dataset and bowtie indices. 

Open the `config.txt` file in the scripts directory. You will need to set the entries enclosed with square brackets [] with files obtained as detailed below. To begin, set `code_directory` to the directory where you unpacked the deFuse code. Set `dataset_directory` to the directory you created for storing the reference dataset and bowtie indices. 

### External Tools

deFuse relies on other publically available tools as part of its pipeline. Some of these tools are not included with the deFuse download. Obtain these tools as detailed below. 

Download samtools: 

The latest version of samtools can be downloaded from the [samtools sourceforge site]](https://sourceforge.net/projects/samtools/files/samtools). 

Set the `samtools_bin` entry in `config.txt` to the fully qualified paths of the `samtools` binary.

Download bowtie: 

The latest version of bowtie can be downloaded from the [bowtie sourceforge site](http://sourceforge.net/projects/bowtie-bio/files/bowtie/). deFuse has been tested on version `0.12.5`. 

Set the `bowtie_bin` and `bowtie_build_bin` entries in `config.txt` to the fully qualified paths of the `bowtie` and `bowtie-build` binaries. 

Download blat and faToTwoBit: 

The latest blat tool suite can be downloaded from the [ucsc website](http://hgdownload.cse.ucsc.edu/admin/exe/). Download `blat` and `faToTwoBit` and set the `blat_bin` and `fatotwobit_bin` entries in `config.txt` to the fully qualified paths of the `blat` and `faToTwoBit` binaries. 

Download GMAP: 

The latest version of GMAP can be downloaded from the [gmap site](http://research-pub.gene.com/gmap/). Build with a default configuration. Do not worry about the `--with-gmapdb` build flag, deFuse will request a specific directory for the database anyway. 

Download R: 

The latest version of R can be downloaded from the [R project website](http://www.r-project.org/). Install R and then locate the R and Rscript executables, and set the `r_bin` and `rscript_bin` entries in `config.txt` to the path of those executables. 

Install the ada package. Run R, then at the prompt type `install.packages("ada")`

### Reference Dataset

The reference dataset setup process has been simplified as of deFuse 0.6.0, and deFuse now automatically downloads all required files. 

The `create_reference_dataset.pl` script will download the genome and other source files, and build any derivative files including bowtie indices, gmap indices, and 2bit files. Run the following command. Expect this step to take at least 12 hours. 

```
create_reference_dataset.pl -c config.txt
```

## Input data formats

deFuse now takes only paired fastq as input. To help create paired fastq files, the scripts `fq_all2std.pl`, `qseq2fastq.pl` are available for export and qseq formats, and the tool `bamfastq` is available for bam files. 

## How to run

Running the deFuse pipeline is simple as running a single script. To run the deFuse pipeline, run `defuse.pl` from the scripts directory with the appropriate command line parameters. The parameters are listed below: 

```
-c, --config
Configuration Filename

-1, --1fastq 
Fastq filename 1

-2, --2fastq 
Fastq filename 2

-o, --output 
Output Directory

-n, --name 
Library Name

-l, --local 
Job Local Directory (default: Output Directory)

-s, --submit 
Submitter Type (default: direct)

-p, --parallel 
Maximum Number of Parallel Jobs (default: 1)
```

In the simplest case, deFuse can be run as follows: 

```
defuse.pl -c config.txt -1 reads1.fq -2 reads2.fq -o output_dir
```

With the above parameters, deFuse will use reads from the files `reads1.fq` and `reads2.fq` the output will be in the `output_dir` directory. 

_**note: the output directory should be different from the directory containing `reads1.fq` and `reads2.fq`**_

The above example will not be the fastest way to run deFuse. Given a machine with multiple processes, 8 for example, run deFuse as follows: 

```
defuse.pl -c config.txt -1 reads1.fq -2 reads2.fq -o output_dir -p 8
```

If you have access to a cluster, you may be able to run deFuse as follows for a sun grid engine (SGE) cluster: 

```
defuse.pl -c config.txt -1 reads1.fq -2 reads2.fq -o output_dir -s sge
```

or as follows for a portable batch system (PBS) cluster: 

```
defuse.pl -c config.txt -1 reads1.fq -2 reads2.fq -o output_dir -s pbs
```

or as follows for a LSF cluster: 

```
defuse.pl -c config.txt -1 reads1.fq -2 reads2.fq -o output_dir -s lsf
```

In many cases it is beneficial to store intermediate results on a local disk rather than a network share. This can be done using the `--local` command line parameters as follows: 

```
defuse.pl -c config.txt -1 reads1.fq -2 reads2.fq -o output_dir -s lsf -l /localdisk
```

to specify that intermediate files be stored on at /localdisk. 

## Output

The output directory specified on the command line of `defuse.pl` will contain the files `results.tsv`, `results.filtered.tsv`, and `results.classify.tsv`. All three files have the same format, though `results.classify.tsv` has a probability column from the application of the classifier to `results.tsv`, and `results.filtered.tsv` has been filtered according to the threshold probability as set in `config.tsv`. The file format is tab delimited with one prediction per line, and the following fields per prediction. 

### Identification

  * **cluster_id** : random identifier assigned to each prediction
  * **library_name** : library name given on the command line of deFuse
  * **gene1** : ensembl id of gene 1
  * **gene2** : ensembl id of gene 2
  * **gene_name1** : name of gene 1
  * **gene_name2** : name of gene 2

### Evidence

  * **break_predict** : breakpoint prediction method, denovo or splitr, that is considered most reliable
  * **concordant_ratio** : proportion of spanning reads considered concordant by blat
  * **denovo_min_count** : minimum kmer count across denovo assembled sequence
  * **denovo_sequence** : fusion sequence predicted by debruijn based denovo sequence assembly
  * **denovo_span_pvalue** : p-value, lower values are evidence the prediction is a false positive
  * **gene_align_strand1** : alignment strand for spanning read alignments to gene 1
  * **gene_align_strand2** : alignment strand for spanning read alignments to gene 2
  * **min_map_count** : minimum of the number of genomic mappings for each spanning read
  * **max_map_count** : maximum of the number of genomic mappings for each spanning read
  * **mean_map_count** : average of the number of genomic mappings for each spanning read
  * **num_multi_map** : number of spanning reads that map to more than one genomic location
  * **span_count** : number of spanning reads supporting the fusion
  * **span_coverage1** : coverage of spanning reads aligned to gene 1 as a proportion of expected coverage
  * **span_coverage2** : coverage of spanning reads aligned to gene 2 as a proportion of expected coverage
  * **span_coverage_min** : minimum of span_coverage1 and span_coverage2
  * **span_coverage_max** : maximum of span_coverage1 and span_coverage2
  * **splitr_count** : number of split reads supporting the prediction
  * **splitr_min_pvalue** : p-value, lower values are evidence the prediction is a false positive
  * **splitr_pos_pvalue** : p-value, lower values are evidence the prediction is a false positive
  * **splitr_sequence** : fusion sequence predicted by split reads
  * **splitr_span_pvalue** : p-value, lower values are evidence the prediction is a false positive

### Annotation

  * **adjacent** : fusion between adjacent genes
  * **altsplice** : fusion likely the product of alternative splicing between adjacent genes
  * **break_adj_entropy1** : di-nucleotide entropy of the 40 nucleotides adjacent to the fusion splice in gene 1
  * **break_adj_entropy2** : di-nucleotide entropy of the 40 nucleotides adjacent to the fusion splice in gene 2
  * **break_adj_entropy_min** : minimum of break_adj_entropy1 and break_adj_entropy2
  * **breakpoint_homology** : number of nucleotides at the fusion splice that align equally well to gene 1 or gene 2
  * **breakseqs_estislands_percident** : maximum percent identity of fusion sequence alignments to est islands
  * **cdna_breakseqs_percident** : maximum percent identity of fusion sequence alignments to cdna
  * **deletion** : fusion produced by a genomic deletion
  * **est_breakseqs_percident** : maximum percent identity of fusion sequence alignments to est
  * **eversion** : fusion produced by a genomic eversion
  * **exonboundaries** : fusion splice at exon boundaries
  * **expression1** : expression of gene 1 as number of concordant pairs aligned to exons
  * **expression2** : expression of gene 2 as number of concordant pairs aligned to exons
  * **gene_chromosome1** : chromosome of gene 1
  * **gene_chromosome2** : chromosome of gene 2
  * **gene_end1** : end position for gene 1
  * **gene_end2** : end position for gene 2
  * **gene_location1** : location of breakpoint in gene 1
  * **gene_location2** : location of breakpoint in gene 2
  * **gene_start1** : start of gene 1
  * **gene_start2** : start of gene 2
  * **gene_strand1** : strand of gene 1
  * **gene_strand2** : strand of gene 2
  * **genome_breakseqs_percident** : maximum percent identity of fusion sequence alignments to genome
  * **genomic_break_pos1** : genomic position in gene 1 of fusion splice / breakpoint
  * **genomic_break_pos2** : genomic position in gene 2 of fusion splice / breakpoint
  * **genomic_strand1** : genomic strand in gene 1 of fusion splice / breakpoint, retained sequence upstream on this strand, breakpoint is downstream
  * **genomic_strand2** : genomic strand in gene 2 of fusion splice / breakpoint, retained sequence upstream on this strand, breakpoint is downstream
  * **interchromosomal** : fusion produced by an interchromosomal translocation
  * **interrupted_index1** : ratio of coverage before and after the fusion splice / breakpoint in gene 1
  * **interrupted_index2** : ratio of coverage before and after the fusion splice / breakpoint in gene 2
  * **inversion** : fusion produced by genomic inversion
  * **orf** : fusion combines genes in a way that preserves a reading frame
  * **probability** : probability produced by classification using adaboost and example positives/negatives (only given in results.classified.tsv)
  * **read_through** : fusion involving adjacent potentially resulting from co-transcription rather than genome rearrangement
  * **repeat_proportion1** : proportion of the spanning reads in gene 1 that span a repeat region
  * **repeat_proportion2** : proportion of the spanning reads in gene 2 that span a repeat region
  * **max_repeat_proportion** : max of repeat_proportion1 and repeat_proportion2
  * **splice_score** : number of nucleotides similar to GTAG at fusion splice
  * **num_splice_variants** : number of potential splice variants for this gene pair
  * **splicing_index1** : number of concordant pairs in gene 1 spanning the fusion splice / breakpoint, divided by number of spanning reads supporting the fusion with gene 2
  * **splicing_index2** : number of concordant pairs in gene 2 spanning the fusion splice / breakpoint, divided by number of spanning reads supporting the fusion with gene 1

## Supporting Reads

Supporting reads for a given fusion prediction can be obtained using the `get_reads.pl` script. For a prediction with `cluster_id=123` run the following command to obtain the spanning and split reads supporting the fusion. 

```
get_reads.pl -c config.txt -o output_dir -i 123
```

Split reads are output first as demarcated by the heading `Split Reads:`. The third line of the output corresponds to the approximate sequence of the fusion boundary to which potential split reads are aligned. Subsequent lines alternate between showing a read id and split alignment. 

Spanning reads are output second as demarcated by the heading `Spanning Reads:`. The spanning read information is tab separated and provides information about the alignments of supporting spanning reads. The fields are: 

  * read id
  * encoded read end (0 -> 1, 1 -> 2)
  * chromosome of alignment
  * strand of alignment
  * start of alignment
  * end of alignment

To obtain the original sequence for spanning reads, search the reads.1.fastq file and reads.2.fastq file in the output directory for the fastq entry with the associated read id. 

### Caveats

You cannot run get_reads.pl if you have removed the temporary files from the output folder. The following files are required: 

  * concordant.read.stats
  * clusters.sc
  * reads.split.catalog
  * job/reads.*.improper.sam
  * job/reads.*.spanning.filelist
  * job/reads.*.spanning/*
