#! /usr/bin/env bash

mkdir -p report_data

# Copy Isosceles files
cp -f isosceles_results/se_transcript_strict.rds \
  report_data/isosceles_se_transcript.rds
for sample_id in LIB5432315_SAM24385458 LIB5432316_SAM24385459 LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
cp -f isosceles_results/${sample_id}_se_transcript.rds \
  report_data/isosceles_${sample_id}_se_transcript.rds
done

# Copy IsoQuant files
for sample_id in LIB5432309_SAM24385452 LIB5432310_SAM24385453 \
LIB5432311_SAM24385454 LIB5432312_SAM24385455 LIB5432313_SAM24385456 \
LIB5432314_SAM24385457 LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
cp -f isoquant_results/${sample_id}/00_${sample_id}.fastq/00_${sample_id}.fastq.transcript_tpm.tsv \
  report_data/isoquant_${sample_id}_transcript_tpm.tsv
done

# Copy FLAMES files
for sample_id in LIB5432309_SAM24385452 LIB5432310_SAM24385453 \
LIB5432311_SAM24385454 LIB5432312_SAM24385455 LIB5432313_SAM24385456 \
LIB5432314_SAM24385457 LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
cp -f flames_results/${sample_id}/transcript_count.csv.gz \
  report_data/flames_${sample_id}_transcript_count.csv.gz
done

# Copy Sicelore files
for sample_id in LIB5432309_SAM24385452 LIB5432310_SAM24385453 \
LIB5432311_SAM24385454 LIB5432312_SAM24385455 LIB5432313_SAM24385456 \
LIB5432314_SAM24385457 LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
cp -f sicelore_results/${sample_id}/sicelore_isomatrix.txt \
  report_data/sicelore_${sample_id}_sicelore_isomatrix.txt
done

# Copy bambu files
for sample_id in LIB5432315_SAM24385458 LIB5432316_SAM24385459 LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
cp -f bambu_results/${sample_id}_counts_transcript.txt \
  report_data/bambu_${sample_id}_counts_transcript.txt
done

# Copy Flair files
for sample_id in LIB5432315_SAM24385458 LIB5432316_SAM24385459 LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
cp -f flair_results/${sample_id}/flair.counts_matrix.tsv.counts.tsv \
  report_data/flair_${sample_id}_counts_matrix.tsv
done

# Copy NanoCount files
for sample_id in LIB5432315_SAM24385458 LIB5432316_SAM24385459 LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
cp -f nanocount_results/${sample_id}.tsv \
  report_data/nanocount_${sample_id}.tsv
done

# Copy LIQA files
for sample_id in LIB5432315_SAM24385458 LIB5432316_SAM24385459 LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
cp -f liqa_results/${sample_id}.tsv \
  report_data/liqa_${sample_id}.tsv
done

# Copy ESPRESSO files
for sample_id in LIB5432315_SAM24385458 LIB5432316_SAM24385459 LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
cp -f espresso_results/${sample_id}/samples_N2_R0_abundance.esp \
  report_data/espresso_${sample_id}_abundance.esp
done
