#! /usr/bin/env bash

mkdir -p report_data

# Copy Sequin reference data files
cp -f reference_data/sequin_transcripts.tsv \
  report_data/sequin_transcripts.tsv

# Copy Isosceles files
for sample_id in mixA_ont mixB_ont; do
cp -f isosceles_results/${sample_id}_se_transcript.rds \
  report_data/isosceles_${sample_id}_se_transcript.rds
done

# Copy IsoQuant files
for sample_id in mixA_ont mixB_ont; do
cp -f isoquant_results/${sample_id}/00_${sample_id}.fastq/00_${sample_id}.fastq.transcript_tpm.tsv \
  report_data/isoquant_${sample_id}_transcript_tpm.tsv
done

# Copy bambu files
for sample_id in mixA_ont mixB_ont; do
cp -f bambu_results/${sample_id}_counts_transcript.txt \
  report_data/bambu_${sample_id}_counts_transcript.txt
done

# Copy Flair files
for sample_id in mixA_ont mixB_ont; do
cp -f flair_results/${sample_id}/flair.counts_matrix.tsv.counts.tsv \
  report_data/flair_${sample_id}_counts_matrix.tsv
done

# Copy NanoCount files
for sample_id in mixA_ont mixB_ont; do
cp -f nanocount_results/${sample_id}.tsv \
  report_data/nanocount_${sample_id}.tsv
done

# Copy LIQA files
for sample_id in mixA_ont mixB_ont; do
cp -f liqa_results/${sample_id}.tsv \
  report_data/liqa_${sample_id}.tsv
done

# Copy ESPRESSO files
for sample_id in mixA_ont mixB_ont; do
cp -f espresso_results/${sample_id}/samples_N2_R0_abundance.esp \
  report_data/espresso_${sample_id}_abundance.esp
done
