#! /usr/bin/env bash

# Set the number of CPUs/threads for the analysis
export ncpu=1

# Set the sample ID
export sample_id="cdna"

# Download SIRV annotation data
mkdir -p reference_data
wget -O reference_data/sirvome.fasta \
  http://s3.amazonaws.com/nanopore-human-wgs/rna/referenceFastaFiles/sirv/SIRVome_isoforms_ERCCs_170612a.fasta
wget -O reference_data/sirvome_correct_annotations.gtf \
  http://s3.amazonaws.com/nanopore-human-wgs/rna/referenceFastaFiles/sirv/SIRVome_isoforms_ERCCs_Lot001485_C_170612a.gtf
wget -O reference_data/sirvome_insufficient_annotations.gtf \
  http://s3.amazonaws.com/nanopore-human-wgs/rna/referenceFastaFiles/sirv/SIRVome_isoforms_ERCCs_Lot001485_I_170612a.gtf
wget -O reference_data/sirvome_over_annotations.gtf \
  http://s3.amazonaws.com/nanopore-human-wgs/rna/referenceFastaFiles/sirv/SIRVome_isoforms_ERCCs_Lot001485_O_170612a.gtf

# Download the SIRV BAM file
mkdir -p bam
wget -O bam/sirvome_cdna.bam \
  https://s3.amazonaws.com/nanopore-human-wgs/rna/bamFiles/NA12878-cDNA-1D.pass.dedup.fastq.SIRVome.minimap2.sorted.bam
wget -O bam/sirvome_cdna.bam.bai \
  https://s3.amazonaws.com/nanopore-human-wgs/rna/bamFiles/NA12878-cDNA-1D.pass.dedup.fastq.SIRVome.minimap2.sorted.bam.bai

# Run Isosceles
singularity exec ../singularity/isosceles.sif Rscript run_isosceles.R

# Run bambu
singularity exec ../singularity/bambu.sif Rscript run_bambu.R

# Run IsoQuant
conda activate isosceles_isoquant
mkdir -p isoquant_results
for annotations_id in insufficient over; do
isoquant.py \
  --threads $ncpu \
  --data_type nanopore \
  --reference reference_data/sirvome.fasta \
  --genedb reference_data/sirvome_${annotations_id}_annotations.gtf \
  --bam bam/sirvome_${sample_id}.bam \
  -o isoquant_results/${sample_id}_${annotations_id}
cat reference_data/sirvome_${annotations_id}_annotations.gtf > \
  isoquant_results/${sample_id}_${annotations_id}/extended_annotations.gtf
cat isoquant_results/${sample_id}_${annotations_id}/00_sirvome_${sample_id}/00_sirvome_${sample_id}.transcript_models.gtf \
  | grep "nic" >> \
  isoquant_results/${sample_id}_${annotations_id}/extended_annotations.gtf
isoquant.py \
  --threads $ncpu \
  --no_model_construction \
  --data_type nanopore \
  --reference reference_data/sirvome.fasta \
  --genedb isoquant_results/${sample_id}_${annotations_id}/extended_annotations.gtf \
  --bam bam/sirvome_${sample_id}.bam \
  -o isoquant_results/${sample_id}_${annotations_id}_merged
done
conda deactivate
