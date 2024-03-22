#! /usr/bin/env bash

mkdir -p report_data

# Set the sample ID
export sample_id="cdna"

# Copy reference annotations data parsed by Isosceles
for annotations_id in correct insufficient over; do
cp -f isosceles_results/reference_data/annotations_${annotations_id}.rds \
  report_data/annotations_${annotations_id}.rds
done

# Copy Isosceles files
for annotations_id in insufficient over; do
cp -f isosceles_results/${sample_id}_${annotations_id}_se_transcript.rds \
  report_data/isosceles_${sample_id}_${annotations_id}_se_transcript.rds
done

# Copy IsoQuant files
for annotations_id in insufficient over; do
cp -f isoquant_results/${sample_id}_${annotations_id}_merged/00_sirvome_${sample_id}/00_sirvome_${sample_id}.transcript_tpm.tsv \
  report_data/isoquant_${sample_id}_${annotations_id}_transcript_tpm.tsv
cp -f isoquant_results/${sample_id}_${annotations_id}/extended_annotations.gtf \
  report_data/isoquant_${sample_id}_${annotations_id}_extended_annotations.gtf
done

# Copy bambu files
for annotations_id in insufficient over; do
cp -f bambu_results/${sample_id}_${annotations_id}_counts_transcript.txt \
  report_data/bambu_${sample_id}_${annotations_id}_counts_transcript.txt
cp -f bambu_results/${sample_id}_${annotations_id}_extended_annotations.gtf \
  report_data/bambu_${sample_id}_${annotations_id}_extended_annotations.gtf
done
