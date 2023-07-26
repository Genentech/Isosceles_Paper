#! /usr/bin/env bash

mkdir -p report_data

# Set the sample ID
export sample_id="truncated_scrnaseq"

# Copy Isosceles files
cp -f isosceles_results/se_transcript.rds \
  report_data/isosceles_se_transcript.rds
cp -f isosceles_results/se_pseudobulk_transcript.rds \
  report_data/isosceles_se_pseudobulk_transcript.rds

# Copy IsoQuant files
cp -f isoquant_results/00_${sample_id}/00_${sample_id}.transcript_grouped_tpm.tsv \
  report_data/isoquant_transcript_grouped_tpm.tsv

# Copy FLAMES files
cp -f flames_results/transcript_count.csv.gz \
  report_data/flames_transcript_count.csv.gz

# Copy Sicelore files
cp -f sicelore_results/sicelore_isomatrix.txt \
  report_data/sicelore_isomatrix.txt
