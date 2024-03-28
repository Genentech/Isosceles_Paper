#! /usr/bin/env bash

# Set the number of CPUs/threads for the analysis
export ncpu=1

# Download ENCODE Illumina transcript quantification files
mkdir -p data
wget -O data/ENCFF485OUK.tsv \
  https://www.encodeproject.org/files/ENCFF485OUK/@@download/ENCFF485OUK.tsv

# Download the Pacbio FASTQ file
mkdir -p fastq
wget -O fastq/ENCFF450VAU.fastq.gz \
  https://www.encodeproject.org/files/ENCFF450VAU/@@download/ENCFF450VAU.fastq.gz

# Download the Nanopore FASTQ file
mkdir -p fastq
wget -O fastq/NA12878_ONT.fastq \
  https://s3.amazonaws.com/nanopore-human-wgs/rna/fastq/NA12878-cDNA-1D.pass.dedup.fastq
gzip fastq/NA12878_ONT.fastq

# Align Pacbio reads to the reference genome with minimap2
mkdir -p bam
../software/bin/minimap2 -t $ncpu \
  -ax splice:hq --secondary=no \
  --junc-bed ../reference_data/known_introns.bed --junc-bonus 15 \
  ../reference_data/genome.fasta \
  fastq/ENCFF450VAU.fastq.gz \
  | samtools sort -o bam/ENCFF450VAU.bam
samtools index bam/ENCFF450VAU.bam

# Align Nanopore reads to the reference genome with minimap2
mkdir -p bam
../software/bin/minimap2 -t $ncpu \
  -ax splice --secondary=no \
  --junc-bed ../reference_data/known_introns.bed --junc-bonus 15 \
  ../reference_data/genome.fasta \
  fastq/NA12878_ONT.fastq.gz \
  | samtools sort -o bam/NA12878_ONT.bam
samtools index bam/NA12878_ONT.bam

# Run Isosceles
singularity exec ../singularity/isosceles.sif Rscript run_isosceles.R

# Run bambu
singularity exec ../singularity/bambu.sif Rscript run_bambu.R

# Run IsoQuant (Pacbio data)
conda activate isosceles_isoquant
mkdir -p isoquant_results
isoquant.py \
  --threads $ncpu \
  --no_model_construction \
  --data_type pacbio_ccs \
  --reference ../reference_data/genome.fasta \
  --genedb ../reference_data/ensembl_90.gtf \
  --fastq fastq/ENCFF450VAU.fastq.gz \
  -o isoquant_results/ENCFF450VAU
conda deactivate

# Run IsoQuant (Nanopore data)
conda activate isosceles_isoquant
mkdir -p isoquant_results
isoquant.py \
  --threads $ncpu \
  --no_model_construction \
  --data_type nanopore \
  --reference ../reference_data/genome.fasta \
  --genedb ../reference_data/ensembl_90.gtf \
  --fastq fastq/NA12878_ONT.fastq.gz \
  -o isoquant_results/NA12878_ONT
conda deactivate
