# TrAEL-seq

A collection of scripts for and information about TrAEL-seq. Last update 21 July 2020

### TrAEL-seq pre-processing (UMI and Poly-T handling)

Raw TrAEL-seq FastQ reads are expected to have the following structure:

```
barcode (UMI) (8bp)    //    PolyT     //      Insert
```

**Step 1:**

The script `TrAELseq_preprocessing.py` removes the first 8bp (UMI) of a read and adds the sequence to the end of the readID (separated with a colon, like so:`:UMISEQUENCE`). The quality information is discarded. Empty spaces in the readID are replaced with `_` to preserve the UMI sequence after mapping.

**Step 2:**

After moving the UMI sequences, the script looks for up to 3 T at the start of the sequence, and removes those. Sequences with more than 3 Ts at the 5' end are clipped a maximum of 3 TTT. This pre-processing script requires Python 3.

In its current form the script takes all FastQ files in a folder, and applies UMI as well as Poly-T handling automatically (sequentially). Usage is:

```
./TrAELseq_preprocessing.py
```
### Adapter-/quality trimming

Following pre-processing, reads need to undergo adapter- and quality as per usual. A simple [Trim Galore](https://github.com/FelixKrueger/TrimGalore) run like this should do the trick:

```
trim_galore file_UMIed.fastq.gz
```

### Alignment

UMI-pre-processed and adapter trimmed files were then aligned to the respective genome using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) using local alignments (option: `--local`).

### Deduplication

Finally, alignment results files [BAM] were then de-duplicated using [UmiBam](https://github.com/FelixKrueger/Umi-Grinder). This takes the mapping position as well as the UMI sequence into account (as perfect matches only).

# Data analysis

### TrAEL-seq_ReadCountPlots.R - 

This is a basic script for displaying read count quantification from TrAEL-seq over specified regions of the genome. As input it uses an annotated probe report exported from Seqmonk following quantification of reads (truncated to 1bp) using running windows of any size. Running this whole script should output a pdf file for a specified sample.

### TrAEL-seq_RFDPlots.R - 

This is a basic script for calculating and displaying RFD (replication fork directionality) from TrAEL-seq over specified regions of the genome. As input it uses annotated probe reports from Seqmonk containing forward and reverse read counts (truncated to 1bp) for running windows of any size. Two types of plot can be output (separated by dashed lines in the script) - red and blue dot plots for individual samples, and line plots for multiple samples. 



### eccDNA-seq pre-processing (UMI and single-T handling)

A new eccDNA-seq protocol attempting to introduce UMIs produces reads with the following structure:

```
barcode (UMI) (8bp)    //    T     //      Insert
```

Processing by `eccDNAseq_preprocessing.py` removed the UMI and write the sequence in the readIDs, and removes the next base which should be a T in 100% of cases.
