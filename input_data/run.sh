#! /usr/bin/env bash

# Download MinION bulk RNA-Seq FASTQ files from NCBI GEO (Dataset ID: TBD),
# then save them in the fastq_ont/ directory.
# Following files need to be created:
# * fastq_ont/LIB5427896_SAM24376275.fastq.gz
# * fastq_ont/LIB5427897_SAM24376276.fastq.gz

# Download Promethion bulk RNA-Seq FASTQ files from NCBI GEO (Dataset ID: TBD),
# then save them in the fastq_ont/ directory.
# Following files need to be created:
# * fastq_ont/LIB5432309_SAM24385452.fastq.gz
# * fastq_ont/LIB5432310_SAM24385453.fastq.gz
# * fastq_ont/LIB5432311_SAM24385454.fastq.gz
# * fastq_ont/LIB5432312_SAM24385455.fastq.gz
# * fastq_ont/LIB5432313_SAM24385456.fastq.gz
# * fastq_ont/LIB5432314_SAM24385457.fastq.gz
# * fastq_ont/LIB5432315_SAM24385458.fastq.gz
# * fastq_ont/LIB5432316_SAM24385459.fastq.gz

# Download Nanopore scRNA-Seq FASTQ files from NCBI GEO (Dataset ID: TBD),
# then save them in the fastq_ont/ directory.
# Following files need to be created:
# * fastq_ont/LIB5445493_SAM24404003.fastq.gz

# Download Illumina scRNA-Seq BAM files from NCBI GEO (Dataset ID: TBD),
# then save them in the cellranger/ directory.
# Following files need to be created:
# * cellranger/possorted_genome_bam.bam
# Build the BAM file index using the following command:
# samtools index cellranger/possorted_genome_bam.bam

# Download simulated Nanopore FASTQ files
wget https://zenodo.org/record/8180696/files/truncated_bulk_rnaseq.fastq.gz
wget https://zenodo.org/record/8180696/files/truncated_scrnaseq.fastq.gz

# Download the simulated Nanopore scRNA-Seq BAM file
wget https://zenodo.org/record/8180696/files/truncated_scrnaseq.bam
samtools index bam_sim/truncated_scrnaseq.bam
