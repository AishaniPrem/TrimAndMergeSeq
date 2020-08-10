# TrimAndMergeSeq
This program is used to Trim off the primers of the provided demultiplexed Illumina paired-end reads and merge the reads. This program uses the DADA2 pakage of R to trim 16S sequences (V4-V5 region) with primers size 27.

## To run the code.
On the command line run the following scripts.

`sudo bash install.txt`

`Rscript TrimAndMerge.r path-to-folder-contning-demultiplexed-fastqfiles`

## Output files.
* Folder "Merged" containing the merged and trimmed files
* "SequenceAbundanceSummary.csv" file contains sequence abundance distribution for each sequence in each sample.
* "SequenceID.csv" file contains the list of sequences detected along with the assigned sequence id.
