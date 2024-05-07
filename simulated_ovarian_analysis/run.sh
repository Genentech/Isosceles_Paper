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

# Run IsoQuant (scRNA-Seq, Rep1)
conda activate isosceles_isoquant
mkdir -p isoquant_results
for sample_id in SK-OV-3 IGROV-1 OVMANA OVKATE OVTOKO COV362; do
isoquant.py \
  --threads $ncpu \
  --no_model_construction \
  --data_type nanopore \
  --reference ../reference_data/genome.fasta \
  --genedb ../reference_data/simulated_transcript_annotations.gtf \
  --bam ../input_data/bam_sim/${sample_id}_Rep1.bam \
  --read_group tag:BC \
  -o isoquant_results/${sample_id}_Rep1
done
conda deactivate

# Run IsoQuant (bulk RNA-Seq, Rep2)
conda activate isosceles_isoquant
mkdir -p isoquant_results
for sample_id in SK-OV-3 IGROV-1 OVMANA OVKATE OVTOKO COV362; do
isoquant.py \
  --threads $ncpu \
  --no_model_construction \
  --data_type nanopore \
  --reference ../reference_data/genome.fasta \
  --genedb ../reference_data/simulated_transcript_annotations.gtf \
  --fastq ../input_data/fastq_sim/${sample_id}_Rep2.fastq.gz \
  -o isoquant_results/${sample_id}_Rep2
done
conda deactivate

# Prepare a FASTQ files from simulated BAM files for FLAMES (scRNA-Seq, Rep1)
mkdir -p fastq
for sample_id in SK-OV-3 IGROV-1 OVMANA OVKATE OVTOKO COV362; do
bedtools bamtofastq -i ../input_data/bam_sim/${sample_id}_Rep1.bam \
  -fq fastq/${sample_id}_Rep1_flames.fastq
gzip fastq/${sample_id}_Rep1_flames.fastq
done

# Run FLAMES (scRNA-Seq, Rep1)
conda activate isosceles_flames_deps
mkdir -p flames_results
for sample_id in SK-OV-3 IGROV-1 OVMANA OVKATE OVTOKO COV362; do
../software/git/FLAMES/python/sc_long_pipeline.py \
  -a ../reference_data/simulated_transcript_annotations.gtf \
  -i fastq/${sample_id}_Rep1_flames.fastq.gz \
  --outdir flames_results/${sample_id}_Rep1 \
  --genomefa ../reference_data/genome.fasta \
  --minimap2_dir $project_dir/conda_flames_deps/bin
done
conda deactivate

# Run FLAMES (bulk RNA-Seq, Rep2)
conda activate isosceles_flames_deps
mkdir -p flames_results
for sample_id in SK-OV-3 IGROV-1 OVMANA OVKATE OVTOKO COV362; do
mkdir -p fastq/${sample_id}_Rep2
ln -s ../input_data/fastq_sim/${sample_id}_Rep2.fastq.gz \
  fastq/${sample_id}_Rep2/${sample_id}_Rep2.fastq.gz
../software/git/FLAMES/python/bulk_long_pipeline.py \
  --config_file ../software/configs/config_sclr_nanopore_bulk.json \
  -a ../reference_data/simulated_transcript_annotations.gtf \
  -i fastq/${sample_id}_Rep2 \
  --outdir flames_results/${sample_id}_Rep2 \
  --genomefa ../reference_data/genome.fasta \
  --minimap2_dir $(which minimap2 | xargs dirname)
rm -rf fastq/${sample_id}_Rep2
done
conda deactivate

# Prepare cell barcodes for Sicelore (bulk RNA-Seq, Rep2)
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

# Run Sicelore (scRNA-Seq, Rep1)
mkdir -p sicelore_results
for sample_id in SK-OV-3 IGROV-1 OVMANA OVKATE OVTOKO COV362; do
mkdir -p sicelore_results/${sample_id}_Rep1
java $java_xmx_flag -jar ../software/git/sicelore/Jar/Sicelore-2.0.jar \
  IsoformMatrix \
  I=bam/${sample_id}_Rep1_sicelore_GE.bam \
  OUTDIR=sicelore_results/${sample_id}_Rep1 \
  PREFIX=sicelore \
  ISOBAM=true \
  GENETAG=GE UMITAG=U8 CELLTAG=BC \
  REFFLAT=reference_data/simulated_transcript_annotations.refFlat.txt \
  CSV=../reference_data/cell_barcodes.txt \
  DELTA=2 MAXCLIP=150 METHOD=STRICT AMBIGUOUS_ASSIGN=false
done

# Run Sicelore (bulk RNA-Seq, Rep2)
mkdir -p sicelore_results
for sample_id in SK-OV-3 IGROV-1 OVMANA OVKATE OVTOKO COV362; do
mkdir -p sicelore_results/${sample_id}_Rep2
java $java_xmx_flag -jar ../software/git/sicelore/Jar/Sicelore-2.0.jar \
  IsoformMatrix \
  I=bam/${sample_id}_Rep2_sicelore_GE.bam \
  OUTDIR=sicelore_results/${sample_id}_Rep2 \
  PREFIX=sicelore \
  ISOBAM=true \
  GENETAG=GE UMITAG=U8 CELLTAG=BC \
  REFFLAT=reference_data/simulated_transcript_annotations.refFlat.txt \
  CSV=reference_data/barcodes.csv \
  DELTA=2 MAXCLIP=150 METHOD=STRICT AMBIGUOUS_ASSIGN=false
done
