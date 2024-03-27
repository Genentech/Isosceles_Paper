#! /usr/bin/env bash

# Set the number of CPUs/threads for the analysis
export ncpu=1

# Set the Java max memory flag
export java_xmx_flag="-Xmx1g"

# Align Nanopore reads to the reference genome with minimap2 (bulk RNA-Seq, Rep2)
mkdir -p bam
for sample_id in SK-OV-3 IGROV-1 OVMANA OVKATE OVTOKO COV362; do
../software/bin/minimap2 -t $ncpu \
  -ax splice --secondary=no \
  --junc-bed ../reference_data/simulated_known_introns.bed --junc-bonus 15 \
  ../reference_data/genome.fasta \
  ../input_data/fastq_sim/${sample_id}_Rep2.fastq.gz \
  | samtools sort -o bam/${sample_id}_Rep2.bam
samtools index bam/${sample_id}_Rep2.bam
done

# Run Isosceles
singularity exec ../singularity/isosceles.sif Rscript run_isosceles_bulk.R
singularity exec ../singularity/isosceles.sif Rscript run_isosceles_sc.R

# Prepare a FASTQ files from simulated BAM files for FLAMES (scRNA-Seq, Rep1)
mkdir -p fastq
for sample_id in SK-OV-3 IGROV-1 OVMANA OVKATE OVTOKO COV362; do
bedtools bamtofastq -i ../input_data/bam_sim/${sample_id}_Rep1.bam \
  -fq fastq/${sample_id}_Rep1_flames.fastq
gzip fastq/${sample_id}_Rep1_flames.fastq
done

# Prepare cell barcodes for Sicelore
mkdir -p reference_data
echo AAACCCAAGTATTGCC > reference_data/barcodes.csv

# Prepare the reference annotations (refFlat) file for Sicelore
mkdir -p reference_data
gtfToGenePred -genePredExt -ignoreGroupsWithoutExons \
  ../reference_data/simulated_transcript_annotations.gtf \
  /dev/stdout | \
  awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > \
  reference_data/simulated_transcript_annotations.refFlat.txt

# Tag the BAM files with gene names for Sicelore (scRNA-Seq, Rep1)
for sample_id in SK-OV-3 IGROV-1 OVMANA OVKATE OVTOKO COV362; do
java $java_xmx_flag -jar ../software/git/sicelore/Jar/Sicelore-2.0.jar \
  AddGeneNameTag \
  I=../input_data/bam_sim/${sample_id}_Rep1.bam \
  O=bam/${sample_id}_Rep1_sicelore_GE.bam \
  REFFLAT=reference_data/simulated_transcript_annotations.refFlat.txt
samtools index bam/${sample_id}_Rep1_sicelore_GE.bam
done

# Add cell barcode and UMI tags to the BAM files for Sicelore (bulk RNA-Seq, Rep2)
conda activate isosceles_nanocount
python add_bam_tags.py
conda deactivate

# Tag the BAM files with gene names for Sicelore (bulk RNA-Seq, Rep2)
for sample_id in SK-OV-3 IGROV-1 OVMANA OVKATE OVTOKO COV362; do
java $java_xmx_flag -jar ../software/git/sicelore/Jar/Sicelore-2.0.jar \
  AddGeneNameTag \
  I=bam/${sample_id}_Rep2_sicelore.bam \
  O=bam/${sample_id}_Rep2_sicelore_GE.bam \
  REFFLAT=reference_data/simulated_transcript_annotations.refFlat.txt
samtools index bam/${sample_id}_Rep2_sicelore_GE.bam
done
