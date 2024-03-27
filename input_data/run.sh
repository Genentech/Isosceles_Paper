#! /usr/bin/env bash

# Download Promethion bulk RNA-Seq FASTQ files from NCBI (GEO SuperSeries GSE248118)
fastq-dump --gzip SRR26865809
mv SRR26865809.fastq.gz fastq_ont/LIB5432309_SAM24385452.fastq.gz
fastq-dump --gzip SRR26865808
mv SRR26865808.fastq.gz fastq_ont/LIB5432310_SAM24385453.fastq.gz
fastq-dump --gzip SRR26865807
mv SRR26865807.fastq.gz fastq_ont/LIB5432311_SAM24385454.fastq.gz
fastq-dump --gzip SRR26865806
mv SRR26865806.fastq.gz fastq_ont/LIB5432312_SAM24385455.fastq.gz
fastq-dump --gzip SRR26865805
mv SRR26865805.fastq.gz fastq_ont/LIB5432313_SAM24385456.fastq.gz
fastq-dump --gzip SRR26865804
mv SRR26865804.fastq.gz fastq_ont/LIB5432314_SAM24385457.fastq.gz
fastq-dump --gzip SRR26865803
mv SRR26865803.fastq.gz fastq_ont/LIB5432315_SAM24385458.fastq.gz
fastq-dump --gzip SRR26865802
mv SRR26865802.fastq.gz fastq_ont/LIB5432316_SAM24385459.fastq.gz

# Download MinION bulk RNA-Seq FASTQ files from NCBI (GEO SuperSeries GSE248118)
fastq-dump --gzip SRR26865801
mv SRR26865801.fastq.gz fastq_ont/LIB5427896_SAM24376275.fastq.gz
fastq-dump --gzip SRR26865800
mv SRR26865800.fastq.gz fastq_ont/LIB5427897_SAM24376276.fastq.gz

# Download Nanopore scRNA-Seq FASTQ file from NCBI (GEO SuperSeries GSE248118)
fastq-dump --gzip SRR26865982
mv SRR26865982.fastq.gz fastq_ont/LIB5445493_SAM24404003.fastq.gz

# Download Illumina scRNA-Seq BAM file from NCBI (GEO SuperSeries GSE248118)
wget -O cellranger/possorted_genome_bam.bam \
  https://sra-pub-src-1.s3.amazonaws.com/SRR26865983/LIB5445371_SAM24404003.bam.1
samtools index cellranger/possorted_genome_bam.bam

# Download simulated Nanopore FASTQ files
cd fastq_sim
wget https://zenodo.org/record/8180696/files/truncated_bulk_rnaseq.fastq.gz
wget https://zenodo.org/record/8180696/files/truncated_scrnaseq.fastq.gz
wget https://zenodo.org/record/0000000/files/COV362_Rep2.fastq.gz
wget https://zenodo.org/record/0000000/files/IGROV-1_Rep2.fastq.gz
wget https://zenodo.org/record/0000000/files/OVKATE_Rep2.fastq.gz
wget https://zenodo.org/record/0000000/files/OVMANA_Rep2.fastq.gz
wget https://zenodo.org/record/0000000/files/OVTOKO_Rep2.fastq.gz
wget https://zenodo.org/record/0000000/files/SK-OV-3_Rep2.fastq.gz
cd ..

# Download simulated Nanopore scRNA-Seq BAM files
cd bam_sim
wget https://zenodo.org/record/8180696/files/truncated_scrnaseq.bam
samtools index truncated_scrnaseq.bam
wget https://zenodo.org/record/0000000/files/COV362_Rep1.bam
samtools index COV362_Rep1.bam
wget https://zenodo.org/record/0000000/files/IGROV-1_Rep1.bam
samtools index IGROV-1.bam
wget https://zenodo.org/record/0000000/files/OVKATE_Rep1.bam
samtools index OVKATE_Rep1.bam
wget https://zenodo.org/record/0000000/files/OVMANA_Rep1.bam
samtools index OVMANA_Rep1.bam
wget https://zenodo.org/record/0000000/files/OVTOKO_Rep1.bam
samtools index OVTOKO_Rep1.bam
wget https://zenodo.org/record/0000000/files/SK-OV-3_Rep1.bam
samtools index SK-OV-3_Rep1.bam
cd ..
