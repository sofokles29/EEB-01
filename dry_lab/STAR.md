# Mapping with STAR
To determine where on the brown bear genome our reads originated from, we will align our reads to the reference genome using STAR (Spliced Transcripts Alignment to a Reference). STAR is an aligner designed to specifically address many of the challenges of RNA-seq data mapping using a strategy to account for spliced alignments. We can say that STAR splits each read into chunks, finds where those chunks match the genome, trusts the ones that map uniquely, groups nearby matches, and then connects them into a full read alignment — scoring each possibility to pick the best one.

## Purpose

STAR is shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed, but it is memory intensive. The algorithm achieves this highly efficient mapping by performing a two-step process:

**1. Seed searching**
- For every read that STAR aligns, STAR will search for the longest sequence that exactly matches one or more locations on the reference genome. These longest matching sequences are called the Maximal Mappable Prefixes (MMPs).

- The different parts of the read that are mapped separately are called ‘seeds’ -> the first MMP that is mapped to the genome is called seed1.

- STAR will then search again for only the unmapped portion of the read to find the next longest sequence that exactly matches the reference genome, or the next MMP -> seed2.

- This sequential searching of only the unmapped portions of reads underlies the efficiency of the STAR algorithm.

- STAR uses an uncompressed suffix array (SA) to efficiently search for the MMPs, this allows for quick searching against even the largest reference genomes.

- Other slower aligners use algorithms that often search for the entire read sequence before splitting reads and performing iterative rounds of mapping.

- If STAR does not find an exact matching sequence for each part of the read due to mismatches or indels, the previous MMPs will be extended.

- If extension does not give a good alignment, then the poor quality or adapter sequence (or other contaminating sequence) will be soft clipped.


**2. Clustering, stitching, and scoring**
- The separate seeds are stitched together to create a complete read by first clustering the seeds together based on proximity to a set of ‘anchor’ seeds, or seeds that are not multi-mapping.
Anchor seeds map uniquely to one place in the genome

- Then the seeds are stitched together based on the best alignment for the read (scoring based on mismatches, indels, gaps, etc.).


## Function

It performs ultrafast spliced alignments by mapping sequencing reads, including those spanning exon-exon junctions, to the genome. STAR builds a suffix array index of the genome, enabling rapid alignment of millions of reads while accurately detecting splice junctions crucial for transcriptome studies.

**Understanding Genome Indexing**
1. Load and Process the Genome FASTA
STAR reads your input genome FASTA file(s), which contain the DNA sequences of chromosomes/scaffolds.
It combines these sequences into a single internal representation and records chromosome names and lengths.


2. Generate a Suffix Array
STAR builds a suffix array and related data structures that allow for ultra-fast searching of subsequences (e.g., mapping reads to genome locations).
It also builds and index of the genome to allow for quick and memory-efficient lookups.


3. Incorporate Transcript Annotations (Optional, but Recommended)
If you provide a --sjdbGTFfile, STAR parses the GTF to find exon-exon splice junctions.
For each annotated junction, STAR extracts a piece of the genome around the junction (--sjdbOverhang, typically Read Length - 1).
These sequences are added to a splice junction database, improving STAR’s ability to align spliced RNA-seq reads that cross exon boundaries.


4. Prepare Output Index Files
STAR saves a set of binary and text files in the --genomeDir, which include:
Genome — compressed version of the genome sequence.
SA and SAindex — suffix array and its index.
chrName.txt — chromosome names.
sjdbList.out.tab — list of splice junctions (if GTF provided).
sjdbInfo.txt — info about splice junctions.
Other internal files.


These files are later used in the alignment step when you provide --genomeDir.


## Usage

Researchers use STAR after trimming and quality control of RNA-Seq data. The tool is run by specifying the cleaned FASTQ files, the genome index, and output options, producing aligned read files (typically in BAM format) that can be used for quantifying gene expression, identifying novel transcripts, or differential expression analysis in downstream bioinformatics workflows.


**SLURM script:**

```#!/bin/bash

#SBATCH --job-name=STAR 			# Job name
#SBATCH --partition=128x24				# Partition name
#SBATCH --mail-type=ALL               		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=smallinu@ucsc.edu   	# Where to send mail
#SBATCH --time=0-10:00:00 				# Wall clock time limit in Days-Hours:min:seconds
#SBATCH --ntasks=1                    		# Run a single task
#SBATCH --cpus-per-task=12              	# Use 4 threads for fasterq-dump
#SBATCH --output=/hb/groups/sip_eeb_01/sofie/scripts/logs/STAR_%j.out    # Standard output and error log
#SBATCH --error=/hb/groups/sip_eeb_01/sofie/scripts/logs/STAR_%j.err     # Standard output and error log
#SBATCH --mem=40G                    		# Allocate memory for the job.
#SBATCH --array=1-11					# array job
```

**The command of the SLURM script:**

```STAR \
--runThreadN 12 \
--runMode genomeGenerate\
--genomeDir 2_star/indexed_genome \
--genomeFastaFiles Ursus-arctos.fna \
--sjdbGTFfile Ursus-arctos.gtf \
--sjdbOverhang 100
--outFilterMultimapNmax 1 \
--readFilesIn 1_trim/read_1_val.fastq 1_trim/read_2_val.fastq \
--twopassMode Basic \
--outFileNamePrefix analysis/2_star/${sra}_
--outSAMtype BAM SortedByCoordinate
```

*Explanation of lines we typed in:*

--outFilterMultimapNmax
Read alignments will be output only if the read maps is less or equal than this value, otherwise no alignments will be output

--twopassMode Basic
In the 2-pass mapping job, STAR will map the reads twice. In the 1st pass, the novel junctions will be detected and inserted into the genome indices
In the 2nd pass, all reads will be re-mapped using annotated (from the GTF file) and novel (detected in the 1st pass) junctions.




