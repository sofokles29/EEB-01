# Trimming with Trimgalore
## Purpose
To improve the quality of sequencing data by removing adapter sequences and low-quality bases from reads before downstream analysis.

## Function
It acts as a wrapper tool that combines Cutadapt for trimming adapter sequences and low-quality bases with FastQC for post-trimming quality control. This ensures that only high-quality, clean reads are retained for accurate alignment and quantification. 

## Usage
Researchers use Trim Galore on raw sequencing data (FASTQ files) as an initial preprocessing step. The command specifies input files, adapter type (if needed), and quality thresholds, generating cleaned FASTQ files that can then be aligned to a reference genome with higher accuracy and confidence. It also runs FastQC on trimmed outputs and provides a quality score. We then used MultiQC to aggragate QC data to create a report. 

`multiqc`

We used this command:

`trim_galore --paired -q 24 --fastqc --fastqc_args "--noextract --nogroup --outdir 1_trim/fastqc" --stringency 5 --illumina --length 50
-o 1_trim --clip_R1 12 --clip_R2 12 [path/to/read1] [path/to/read2]`

--paired: Indicates input files are paired-end (two reads per fragment). Read 1 (R1) starts from the 5' end of the fragment, while Read 2 (R2) starts from the 3' end.

-q 24: Trims low-quality bases with Phred scores below 24 from both ends of reads.

--fastqc: Runs FastQC automatically on the trimmed reads.

--fastqc_args "--noextract --nogroup --outdir 1_trim/fastqc":

--noextract: Keeps the FastQC results as .zip files (doesn't unzip them).

--nogroup: Disables per-tile aggregation, which can hide poor quality regions in large files.

--outdir 1_trim/fastqc: Stores FastQC reports in the specified folder.

--stringency 5: Requires a minimum of 5 bp overlap between the adapter and the read for trimming to occur.

--illumina: Looks specifically for Illumina adapter sequences.

--length 50: Discards reads shorter than 50 bp after trimming.

-o 1_trim: Outputs all trimmed FASTQ files to the 1_trim directory.

--clip_R1 12: Removes the first 12 bases of Read 1.

--clip_R2 12: Removes the first 12 bases of Read 2.

[path/to/read1] [path/to/read2]: Input raw FASTQ files from NCBI.
