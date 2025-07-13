# Downloading SRA 
**What is the SRA?**
Sequence Read Archive (SRA) is a giant library of raw DNA and RNA sequencing data

## Purpose
To provide a centralized, publicly accessible repository for storing and sharing raw sequencing data, enabling transparency, reproducibility, and the reuse of genomic data for further research and comparative analyses.

## Function
SRA files are containers for raw sequencing reads. They are highly compressed, hence why we can’t just download fastq files directly. To get data from the SRA, we use two tools:  

`prefetch`
- Downloads .sra files (raw compressed format used by NCBI) → these files are not readable yet

`fasterq-dump`
- Converts .sra files into .fastq files, which are text files of sequencing reads
- Each read has a name, DNA sequence, and quality score


## Usage
A SLURM script is a shell script used to submit jobs to the computer cluster and tells the cluster what resources your job needs and what commands to run.
It automates and manages tasks like…
- Requesting CPUs, memory, and time
- Running your script
- Saving output/error
- Running jobs in parallel aka arrays


`#!/bin/bash

#SBATCH --job-name=myjob           # Name of the job
#SBATCH --mail-type=ALL               # Mail events
#SBATCH --mail-user=USCSID@ucsc.edu   # Where to send mail 
#SBATCH --output=output_%j.txt     # Standard output (%j = job ID)
#SBATCH --error=error_%j.txt       # Standard error log
#SBATCH --time=01:00:00            # Time limit (hh:mm:ss)
#SBATCH --partition=128x24      # Partition/queue name
#SBATCH --ntasks=1                 # Number of tasks → run on a single CPU
#SBATCH --mem=250M              	  # Memory per node → 250 megabytes`


**For SRA**

`#!/bin/bash

#SBATCH --job-name=getSRA    			# Job name
#SBATCH --partition=128x24				# Partition name
#SBATCH --mail-type=ALL               		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=UCSC_ID@ucsc.edu   	# Where to send mail
#SBATCH --time=0-05:00:00 				# Wall clock time limit in Days-Hours:min:seconds
#SBATCH --ntasks=1                    		# Run a single task
#SBATCH --cpus-per-task=4                  	# Use 4 threads for fasterq-dump
#SBATCH --output=scripts/logs/fasterq-dump_%j.out    # Standard output and error log
#SBATCH --error=scripts/logs/fasterq-dump_%j.err     # Standard output and error log
#SBATCH --mem=8G                    		# Allocate memory for the job.
#SBATCH --array=1-11					# array job`





## What is a fastq file?
Files we get from the sequencing companies that include DNA sequences and the quality score that associates with them. They come in specific 4-line format which includes a header line, a line of sequence, a +, and then the quality scores.
