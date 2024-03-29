#! /usr/bin/env bash

mkdir -p report_data

# Copy Isosceles files (Sicelore input)
cp -f isosceles_results/se_transcript.rds \
  report_data/isosceles_se_transcript.rds
cp -f isosceles_results/se_gene.rds \
  report_data/isosceles_se_gene.rds
cp -f isosceles_results/se_pseudobulk_transcript.rds \
  report_data/isosceles_se_pseudobulk_transcript.rds
cp -f isosceles_results/se_igrov_pseudobulk_k_transcript.rds \
  report_data/isosceles_se_igrov_pseudobulk_k_transcript.rds

# Copy Isosceles files (wf-single-cell input)
cp -f isosceles_results_wf/se_transcript.rds \
  report_data/isosceles_wf_se_transcript.rds
cp -f isosceles_results_wf/se_pseudobulk_transcript.rds \
  report_data/isosceles_wf_se_pseudobulk_transcript.rds

# Copy IsoQuant files
cp -f isoquant_results/00_molecules_GE_tags/00_molecules_GE_tags.transcript_grouped_counts.tsv \
  report_data/isoquant_transcript_grouped_tpm.tsv

# Copy FLAMES files
cp -f flames_results/transcript_count.csv.gz \
  report_data/flames_transcript_count.csv.gz

# Copy Sicelore files
cp -f sicelore_results/sicelore_isomatrix.txt \
  report_data/sicelore_isomatrix.txt
