# RNAseq_analysis_brown_bear_adipose
**Hibernation in brown bears (Ursus arctos) includes major physiological changes, enabling them to endure prolonged periods of fasting and immobility. During hibernation, the bears undergo changes such as reductions in body temperature, metabolic suppression, and temporary insulin resistance. Yet, bears remarkably evade the detrimental health effects associated with prolonged inactivity in humans. This project compares blood samples from brown bears at different time points to investigate gene expression changes associated with hibernation using RNA-sequencing (RNAseq). Identifying genes with significant changes across seasons offers a unique perspective on systemic adaptations, including changes in the immune system and metabolic regulation. The findings from this project have implications for human medicine, including advancements in treating metabolic disorders, identifying mediators of metabolic homeostasis, and potential therapeutic approaches for muscle-wasting diseases.**

**The main goal of this project is to identify which genes change between active and hibernating seasons in adipose (fat) and explore why those changes happen.**

# Bioinformatics Pipeline
**Quality trimming (Trim Galore) -> Map transcripts to brown bear genome (STAR) -> Quantify gene counts (featureCounts)**

**The experimental workflow used to investigate gene expression changes in brown bears during different physiological states involves a combination of wet lab procedures and bioinformatics analyses. The process begins with the collection of adipose tissue samples from male and female bears in two conditions: active season and hibernation. From these adipose tissue samples, RNA is extracted and prepared for bulk RNA sequencing (RNA-Seq) to analyze the transcriptome comprehensively.**

**Once RNA-Seq data are generated, the analysis pipeline starts with quality trimming using Trim Galore, which removes low-quality bases and adapter sequences to ensure that only high-quality reads are retained for downstream analysis. These cleaned reads are then aligned to the brown bear reference genome using STAR, a fast and accurate aligner that maps transcript sequences to their corresponding genomic locations. Following alignment, featureCounts is used to quantify gene counts, producing a matrix that records the number of reads aligning to each gene for every sample.**

**The next step is differential gene expression analysis using DESeq2, which statistically identifies genes that are significantly upregulated or downregulated when comparing active versus hibernating states. To enable robust exploratory data analysis, a Variance Stabilizing Transformation (VST) is applied, normalizing the data and stabilizing variance across genes regardless of their expression level. This transformed data is then used to create heatmaps, visualizing gene expression patterns across samples, and to perform Principal Component Analysis (PCA), which reduces data dimensionality to uncover clustering of samples and major sources of variation in gene expression profiles.**

**To interpret the biological significance of differentially expressed genes, enrichment analyses are conducted. Gene Ontology (GO) and KEGG pathway enrichment analysis identify overrepresented biological processes and pathways, providing insight into functional changes between physiological states. Additionally, Gene Set Enrichment Analysis (GSEA) is performed to determine whether predefined sets of genes show statistically significant, coordinated differences in expression between active and hibernating bears, without imposing arbitrary cutoffs for differential expression. Over-representation analysis (ORA) is also included to assess whether specific gene sets are present more frequently than expected by chance within the list of differentially expressed genes.**

**Overall, this integrated workflow combines meticulous sample collection and RNA preparation with a comprehensive bioinformatics and statistical analysis pipeline to uncover meaningful transcriptional differences underlying physiological adaptations in brown bear hibernation.**




