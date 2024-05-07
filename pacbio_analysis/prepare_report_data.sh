#! /usr/bin/env bash

mkdir -p report_data

# Copy ENCODE Illumina files
cp -f data/ENCFF485OUK.tsv \
  report_data/illumina_ENCFF485OUK.tsv

# Copy Isosceles files
for sample_id in ENCFF450VAU NA12878_ONT; do
cp -f isosceles_results/${sample_id}_se_transcript.rds \
  report_data/isosceles_${sample_id}_se_transcript.rds
done

# Copy IsoQuant files
for sample_id in ENCFF450VAU NA12878_ONT; do
cp -f isoquant_results/${sample_id}/00_${sample_id}.fastq/00_${sample_id}.fastq.transcript_tpm.tsv \
  report_data/isoquant_${sample_id}_transcript_tpm.tsv
done

# Copy bambu files
for sample_id in ENCFF450VAU NA12878_ONT; do
cp -f bambu_results/${sample_id}_counts_transcript.txt \
  report_data/bambu_${sample_id}_counts_transcript.txt
done
