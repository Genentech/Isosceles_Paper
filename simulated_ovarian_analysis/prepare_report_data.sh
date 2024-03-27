#! /usr/bin/env bash

mkdir -p report_data

# Copy Isosceles files
for sample_id in SK-OV-3 IGROV-1 OVMANA OVKATE OVTOKO COV362; do
cp -f isosceles_results_bulk/${sample_id}_Rep2_se_transcript.rds \
  report_data/isosceles_${sample_id}_Rep2_se_transcript.rds
cp -f isosceles_results_sc/${sample_id}_Rep1_se_transcript.rds \
  report_data/isosceles_${sample_id}_Rep1_se_transcript.rds
cp -f isosceles_results_sc/${sample_id}_Rep1_se_pseudobulk_transcript.rds \
  report_data/isosceles_${sample_id}_Rep1_se_pseudobulk_transcript.rds
done

# Copy IsoQuant files
for sample_id in SK-OV-3 IGROV-1 OVMANA OVKATE OVTOKO COV362; do
cp -f isoquant_results/${sample_id}_Rep2/00_${sample_id}_Rep2.fastq/00_${sample_id}_Rep2.fastq.transcript_tpm.tsv \
  report_data/isoquant_${sample_id}_Rep2_transcript_tpm.tsv
cp -f isoquant_results/${sample_id}_Rep1/00_${sample_id}_Rep1_sicelore/00_${sample_id}_Rep1_sicelore.transcript_grouped_counts.tsv \
  report_data/isoquant_${sample_id}_Rep1_transcript_grouped_counts.tsv
done

# Copy FLAMES files
for sample_id in SK-OV-3 IGROV-1 OVMANA OVKATE OVTOKO COV362; do
cp -f flames_results/${sample_id}_Rep2/transcript_count.csv.gz \
  report_data/flames_${sample_id}_Rep2_transcript_count.csv.gz
cp -f flames_results/${sample_id}_Rep1/transcript_count.csv.gz \
  report_data/flames_${sample_id}_Rep1_transcript_count.csv.gz
done

# Copy Sicelore files
for sample_id in SK-OV-3 IGROV-1 OVMANA OVKATE OVTOKO COV362; do
cp -f sicelore_results/${sample_id}_Rep2/sicelore_isomatrix.txt \
  report_data/sicelore_${sample_id}_Rep2_sicelore_isomatrix.txt
cp -f sicelore_results/${sample_id}_Rep1/sicelore_isomatrix.txt \
  report_data/sicelore_${sample_id}_Rep1_sicelore_isomatrix.txt
done
