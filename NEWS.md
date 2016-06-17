# News

## Version 0.8.0 is now available

Improvements: 

  * Install via conda
  * Removed paths from config
  * Specify output files on the command line (courtesy Sergey Mitrofanov)
  * Date and time logged for each step (courtesy Sergey Mitrofanov)

Bugfixes:
  * Samtools command line fixes (courtesy Dusan Randjelovic)

## Version 0.6.2 is now available

Bugfixes: 

  * Updated to work with newer version of gmap

## Version 0.6.1 is now available

Improvements: 

  * deFuse now takes 2 fastq files as input, replacing the previous archaic input method 
  * Reimplemented `interrupted_indexN` and `splicing_indexN`
  * Optionally compute `interrupted_indexN` and `splicing_indexN` using `calculate_extra_annotations`

## Version 0.6.0 is now available

Improvements: 

  * Automatic downloading of reference genome and other required files 
  * Much improved speed due to use of GMAP for the majority of post-process alignments 
  * Switched back to adaboost as described in the original paper, the svm was occasionally classifying some highly expressed fusions as false 

Bugfixes: 

  * fixed issue with `merge_read_stats.pl`

Other Changes: 

  * no longer filtering garbage reads since it was slow, and generally not useful 
  * `--mm` bowtie option for memory mapped bowtie index was removed from default options, add to `bowtie_options` if in the config if needed 

Comments: 

  * The results for some longer read datasets are greatly improved by trimming a variable number of low quality bases from the ends of each read 

## Version 0.5.0 is now available

Fixes the following bugs: 

  * Removed overzealous warning: `Fastq error, unable to interpret readid`
  * fixed findcandidatereads crashes 
  * Increased parallelism for faster runtime, especially during split read alignment and spanning read clustering 
  * removed obsolete break_predict column from output 

No dataset rebuild required. 

## Version 0.4.3 is now available

Added a simulated test dataset and deFuse output based on default parameters [here](http://sourceforge.net/projects/defuse/files/simulation/)

Fixes the following bugs: 

  * Inconsequential `Error: expected GTAG and GTC to be same size` message muted 
  * Removed dependence on human ensembl identifiers to allow for use of mouse genome/gene models 
  * Warning message when deFuse fails to import any sequence data 
  * Different read lengths for the same pair no longer fails 
  * Workaround for occasional samtools merge crash 
  * Wrong gene annotation bugfix 
  * Warning for small fragment lengths with discord_read_trim recommendation 

In addition, the script `get_reads.pl` has been added. Pass it the config filename, output directory, and `cluster_id` of a given fusion to view the supporting spanning and split reads. 

No dataset rebuild required. 

## Version 0.4.2 is now available

Fixes the est island filter that was previously broken. Adds 50 A's to each cDNA sequence to allow for alignment of poly A tail reads. 

Requires a dataset rebuild using the `create_reference_dataset.pl` script. 

## Version 0.4.1 is now available

Version 0.4.1 fixes the `too many reference sequences` bug of when using version 0.4.0 with the latest ensembl. 

To start 0.4.1 where 0.4.0 left off, delete `*cluster*` from the output directory and run 0.4.1 on that output directory. There is no need to rebuild the reference dataset. 

## Version 0.4.0 is now available

Version 0.4.0 provides several advantages over version 0.3 defuse: 

  * Fusions between a gene and an intergenic region are now detected 
  * local realignment of predicted sequenced is used to reduce false positives 
  * blat alignment for annotation is split over multiple jobs for speed and lower memory use 
  * as with version 0.3.7, version 0.4.0 should avoid a combinatorial explosion in the number of clusters 

From version 0.4.0 onwards, a prepackaged dataset will be replaced by step by step instructions on how to build your own dataset. The 0.3 dataset will not be compatible with version 0.4.0 and onwards. 

## Version 0.3.7 is now available

Version 0.3.7 fixes speed and memory issues that would occur with some datasets. If deFuse was taking a long time and generating a very large clusters.txt file, this update should help. In order to run the new version on a dataset that failed with version 0.3.6 or lower, first remove all files matching *clusters* in the output directory created by the previous version, then restart version 0.3.7 using the old output directory. The new version should pick up where the old version left off. 

Known issues: 

  * blat may fail on low memory machines, to be fixed soon by splitting up the reference for blat

## Version 0.3.6 is now available

Version 0.3.6 introduces new annotations that are leveraged using the adaboost classifier for slightly higher accuracy. The filtered output is now based on the probability produced by the adaboost classifier. deFuse 0.3.6 updates can be used to quickly reannotate a deFuse 0.3.5 analysis, however the repeats.txt file from the newly posted dataset is required. To reannotate: 

`rm output_directory/annotations.txt` `./reannotate.pl -c config.txt -o output_directory`

Known issues: 

  * we are currently working on speed issues that may be a problem for some larger datasets
  * blat may fail on low memory machines, to be fixed soon by splitting up the reference for blat

## Version 0.3.5 is now available

Version 0.3.5 is a minor update that fixes a number of small bugs and fixes `reannotate.pl` so that it can actually be used to classify the results of any defuse 0.3.X run. This functionality was not working for 0.3.4. Simply run the following to annotate any 0.3.X run and obtain the new annotations and adaboost probability score. It is not necessary to run `reannotate.pl` if you have already run version 0.3.5 of `defuse.pl`. 

`rm output_directory/annotations.txt` `./reannotate.pl -c config.txt -o output_directory -p max_threads`

Known issues: 

  * blat may fail on low memory machines, to be fixed soon by splitting up the reference for blat
  * the method is somewhat sensitive to the max_insert_size parameter, set to 3 standard deviations above expected fragment length, but no higher

## Version 0.3.4 is now available

Version 0.3.4 uses an adaboost classifier trained on 60 true positives and 61 false positives to produce a single probability score for each fusion. The R ada package is required. You can take full advantage of the adaboost classifier for results produced using other 0.3.X versions of defuse by simply running the `reannotate.pl` script ****Edit: this does not work, please update to version 0.3.5****. Doing a full rerun using 0.3.4 should not be necessary. Once again, the dataset has not changed. 

Changes: 

  * adaboost classifier to produce a single probability score for each prediction
  * additional features calculated for each fusion, including boundary sequence di-nucleotide entropy and fusion boundary homology
  * bugfix for long read lengths and short fragment lengths

Known issues: 

  * blat may fail on low memory machines, to be fixed soon by splitting up the reference for blat
  * the method is somewhat sensitive to the max_insert_size parameter, set to 3 standard deviations above expected fragment length, but no higher

## Version 0.3.3 is now available

Version 0.3.3 reworks the split alignments so they are much quicker and do not hit the disk as much as for previous versions. If your deFuse runs are taking a long time, its a good idea to upgrade. The dataset has not changed. 

Changes: 

  * reworked split alignments to be quicker and less disk intensive
  * discord_read_trim in config.txt allows better use of libraries with long read lengths and shorter fragment lengths
  * quality type for bowtie is configurable by setting bowtie_quals in config.txt

Known issues: 

  * blat may fail on low memory machines, to be fixed soon by splitting up the reference for blat
  * the method is somewhat sensitive to the max_insert_size parameter, set to 3 standard deviations above expected fragment length, but no higher

## Version 0.3.2 is now available

If you were having problems previously with deFuse creating very large files, consuming large amounts of memory and taking long amounts of time, it could be because version 0.3.1 and 0.3.0 were not properly filtering IG rearrangements. The problems are exacerbated if your RNA-Seq data was produced from a tumour with high amounts of immune infiltration. This issue is now fixed in version 0.3.2. Note that results from 0.3.1 and 0.3.0 will not be wrong, but will include a superset of what you're interested in (assuming you are not interested in predicting IG rearrangements). 

Version 0.3.2 also provides a number of other annotations that may be interesting for prioritizing fusions for validation or further experiments. 

Changes: 

  * fixed exclusion of IG rearrangements, thereby fixing overly high resource usage by deFuse
  * denovo breakpoint assembly is now optional via a setting in config.txt
  * annotation of the expression of each gene
  * annotation of fusions spliced at exon boundarys corrected
  * annotation of an interruption index added
  * annotation of an fusion splice index added
  * annotation of genomic position and strand of fusion splice added

Known issues: 

  * blat may fail on low memory machines, to be fixed soon by splitting up the reference for blat
  * libraries with long read lengths and short fragment lengths (ie 75bp reads and 150 bp fragments) will not work, to be fixed by chopping large reads

## Version 0.3.1 is now available

Changes: 

  * Two unnecessary and potentially large files were being produced, this is no longer the case

## Version 0.3.0 is now available

This version represents a major change to the way split reads are calculated and has produced more validated predicitons than 0.2.0.  Not backward compatible with 0.2.0. 

Changes: 

  * New split read calculation
  * Parallelized split read calculation
  * Output now includes column headers

## Version 0.2.0 is now available

This is the first official release of deFuse. 
