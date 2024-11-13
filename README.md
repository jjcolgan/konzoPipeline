# konzoPipeline
Automated pipeline for the assembly, curation and binning of metagenomic assemblies 

This pipeline is used to qced NGS data and produce assembled bacterial genomes and metagenomes which are curated with taxonomy as well as a variaty of functional databases.
The workflow was used to assemble and curate a 2.4 TB consisting of NGS data taken from stool. All steps are parralized thus greatly reducing compute time and increasing processing efficiency. 

The workflow utilizes the following packages: snakemake, bowtie2, anvi'o, samtools and megahit. All must be downloaded prior to running, and the arugement 'conda' set to the users conda enviroment for each package. 
