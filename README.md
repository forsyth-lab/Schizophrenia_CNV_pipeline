# Schizophrenia_CNV_pipeline

## Description
This pipeline is an R-based command-line utility designed to streamline the process of merging CNV (Copy Number Variation) calls generated by different algorithms. The tool enhances the reliability and comprehensiveness of CNV detection by integrating results from multiple algorithms(QuantiSNP, PennCNV), applying stringent quality control measures, annotation, and summarization processes to identify and characterize rare CNVs with potential clinical significance.

## Pipeline Workflow
1. Quality Control:

    - Filters out CNVs based on PennCNV confidence scores or QuantiSNP Log Bayes Factor scores.

    - Applies minimum SNP length criteria and evaluates the overlap ratio between PennCNV and QuantiSNP results.

    - Excludes CNVs overlapping with centromeres, telomeres, or segmental duplications in reference genome datasets.

2. Rare CNV Extraction:

    - Drops common CNVs based on gnomAD frequency data

3. CNV Annotation:

    - Intersects CNVs with known loci, including SSD28 and broader neurodevelopmental disorders (NDDs).

    - Annotates all rare CNVs for overlap with exons of any protein-coding gene in the RefSeq database.

    - Total number of genes with one or more exons deleted and number of genes deleted within each of the 18 BrainSpan neurodevelopmental gene-sets described above were summed to generate deletion burden scores for each subject. 

    - LOEUF scores for the canonical transcript for each gene were converted into percentiles across genes, with higher percentiles indicating greater LOF-intolerance. 

4. CNV Summary Statistics

    - Create summary stats for each subject, including total deletion or duplication length, total LOEUF score etc for further analysis.

## Getting Started

- Clone this repository using the following git command:

    `git clone git@github.com:forsyth-lab/Schizophrenia_CNV_pipeline.git`

    Alternatively, download the source files from the github website (`https://github.com/forsyth-lab/Schizophrenia_CNV_pipeline`)
    

- Once R and its dependencies have been installed, run:

    `cd scripts`

    `Rscript CNV_preprocessing.R`