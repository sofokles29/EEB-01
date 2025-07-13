# Quantifying Gene Level Counts with featureCounts

## Purpose
It is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations. It can be used to count both RNA-seq and genomic DNA-seq reads.

## Function

featureCounts takes as input SAM/BAM files and an annotation file including chromosomal coordinates of features. 
It outputs numbers of reads assigned to features (or meta-features). It also outputs stat info for the overall summarization results, including number of successfully assigned reads and number of reads that failed to be assigned due to various reasons (these reasons are included in the stat info).

Each entry in the provided annotation file is taken as a feature (e.g. an exon). A meta-feature is the aggregation of a set of features (e.g. a gene). The featureCounts program uses the gene_id attribute available in the GTF format annotation to group features into meta-features, ie. features belonging to the same meta-feature have the same gene identifier.
featureCounts can count reads at either feature level or at meta-feature level. When summarizing reads at meta-feature level, read counts obtained for features included in the same meta-feature will be added up to yield the read count for the corresponding meta-feature.


## Usage

We used featureCounts after aligning reads to the reference genome (e.g. with STAR). By specifying the alignment files (BAM), annotation file, and feature type (e.g. gene or exon), featureCounts outputs a table of read counts that can then be imported into statistical tools like DESeq2 or to identify differentially expressed genes between experimental conditions.

**SLURM script:**

```#!/bin/bash

#SBATCH --job-name=featureCounts  			# Job name
#SBATCH --partition=128x24				# Partition name
#SBATCH --mail-type=ALL               		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=smallinu@ucsc.edu   	# Where to send mail
#SBATCH --time=0-05:00:00 				# Wall clock time limit in Days-Hours:min:seconds
#SBATCH --ntasks=1                    		# Run a single task
#SBATCH --cpus-per-task=8                 	# Use 4 threads for fasterq-dump
#SBATCH --output=/hb/groups/sip_eeb_01/sofie/scripts/logs/featureCounts_%j.out    # Standard output and error log
#SBATCH --error=/hb/groups/sip_eeb_01/sofie/scripts/logs/featureCounts_%j.err     # Standard output and error log
#SBATCH --mem=8G                    		# Allocate memory for the job.

featureCounts -p -a data/genome/GCF_023065955.2_UrsArc2.0_genomic.gtf -F 'GTF' -g gene_id -t exon -T 8 -o analysis/3_featureCounts/sofie_rawCounts.txt analysis/2_star/*.bam
```


