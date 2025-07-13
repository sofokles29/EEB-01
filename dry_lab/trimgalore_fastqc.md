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

## SLURM script to run TrimGalore
```#!/bin/bash

#SBATCH --job-name=trimgalore  			# Job name
#SBATCH --partition=128x24				# Partition name
#SBATCH --mail-type=ALL               		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=smallinu@ucsc.edu   	# Where to send mail
#SBATCH --time=0-05:00:00 				# Wall clock time limit in Days-Hours:min:seconds
#SBATCH --ntasks=1                    		# Run a single task
#SBATCH --cpus-per-task=1                  	# Use 4 threads for fasterq-dump
#SBATCH --output=/hb/groups/sip_eeb_01/sofie/scripts/logs/trimgalore_%j.out    # Standard output and error log
#SBATCH --error=/hb/groups/sip_eeb_01/sofie/scripts/logs/trimgalore_%j.err     # Standard output and error log
#SBATCH --mem=8G                    		# Allocate memory for the job.
#SBATCH --array=1-11					# array job
```
