#! /usr/bin/env bash

# Set the number of CPUs/threads for the analysis
export ncpu=2

# Set the Java max memory flag
export java_xmx_flag="-Xmx1g"

# Prepare the FASTQ files
mkdir -p fastq
## Downsample Promethion data to 5 million reads
for sample_id in LIB5432309_SAM24385452 LIB5432310_SAM24385453 \
LIB5432311_SAM24385454 LIB5432312_SAM24385455 LIB5432313_SAM24385456 \
LIB5432314_SAM24385457 LIB5432315_SAM24385458 LIB5432316_SAM24385459; do
zcat ../input_data/fastq_ont/${sample_id}.fastq.gz | \
  head -n 20000000 | gzip -c > \
  fastq/${sample_id}.fastq.gz
done
## Prepare symlinks to the MinION FASTQ files
for sample_id in LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
ln -s ../input_data/fastq_ont/${sample_id}.fastq.gz \
  fastq/${sample_id}.fastq.gz
done

# Align Nanopore reads to the reference genome with minimap2
mkdir -p bam
for sample_id in LIB5432309_SAM24385452 LIB5432310_SAM24385453 \
LIB5432311_SAM24385454 LIB5432312_SAM24385455 LIB5432313_SAM24385456 \
LIB5432314_SAM24385457 LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
../software/bin/minimap2 -t $ncpu \
  -ax splice --secondary=no \
  --junc-bed ../reference_data/known_introns.bed --junc-bonus 15 \
  ../reference_data/genome.fasta \
  fastq/${sample_id}.fastq.gz \
  | samtools sort -o bam/${sample_id}.bam
samtools index bam/${sample_id}.bam
done

# Process the reference annotations
mkdir -p reference_data
cat ../reference_data/ensembl_90.gtf | awk '$3 == "exon"' > \
  reference_data/ensembl_90_exon.gtf

# Run Isosceles
singularity exec ../singularity/isosceles.sif Rscript run_isosceles.R

# Run bambu
singularity exec ../singularity/bambu.sif Rscript run_bambu.R

# Run IsoQuant
conda activate isosceles_isoquant
mkdir -p isoquant_results
for sample_id in LIB5432309_SAM24385452 LIB5432310_SAM24385453 \
LIB5432311_SAM24385454 LIB5432312_SAM24385455 LIB5432313_SAM24385456 \
LIB5432314_SAM24385457 LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
isoquant.py \
  --threads $ncpu \
  --no_model_construction \
  --data_type nanopore \
  --reference ../reference_data/genome.fasta \
  --genedb ../reference_data/ensembl_90.gtf \
  --fastq fastq/${sample_id}.fastq.gz \
  -o isoquant_results/${sample_id}
done
conda deactivate

# Run FLAMES
conda activate isosceles_flames_deps
mkdir -p flames_results
for sample_id in LIB5432309_SAM24385452 LIB5432310_SAM24385453 \
LIB5432311_SAM24385454 LIB5432312_SAM24385455 LIB5432313_SAM24385456 \
LIB5432314_SAM24385457 LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
mkdir -p fastq/${sample_id}
ln -s fastq/${sample_id}.fastq.gz \
  fastq/${sample_id}/${sample_id}.fastq.gz
../software/git/FLAMES/python/bulk_long_pipeline.py \
  --config_file ../software/configs/config_sclr_nanopore_bulk.json \
  -a reference_data/ensembl_90_exon.gtf \
  -i fastq/${sample_id} \
  --outdir flames_results/${sample_id} \
  --genomefa ../reference_data/genome.fasta \
  --minimap2_dir $(which minimap2 | xargs dirname)
rm -rf fastq/${sample_id}
done
conda deactivate

# Prepare cell barcodes for Sicelore
mkdir -p reference_data
echo AAACCCAAGTATTGCC > reference_data/barcodes.csv

# Prepare the reference annotations (refFlat) file for Sicelore
mkdir -p reference_data
gtfToGenePred -genePredExt -ignoreGroupsWithoutExons \
  ../reference_data/ensembl_90.gtf \
  /dev/stdout | \
  awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > \
  reference_data/ensembl_90.refFlat.txt

# Add cell barcode and UMI tags to the BAM files for Sicelore
conda activate isosceles_nanocount
python add_bam_tags.py
conda deactivate

# Tag the BAM files with gene names for Sicelore
for sample_id in LIB5432309_SAM24385452 LIB5432310_SAM24385453 \
LIB5432311_SAM24385454 LIB5432312_SAM24385455 LIB5432313_SAM24385456 \
LIB5432314_SAM24385457 LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
java $java_xmx_flag -jar ../software/git/sicelore/Jar/Sicelore-2.0.jar \
  AddGeneNameTag \
  I=bam/${sample_id}_sicelore.bam \
  O=bam/${sample_id}_sicelore_GE.bam \
  REFFLAT=reference_data/ensembl_90.refFlat.txt
samtools index bam/${sample_id}_sicelore_GE.bam
done

# Run Sicelore
for sample_id in LIB5432309_SAM24385452 LIB5432310_SAM24385453 \
LIB5432311_SAM24385454 LIB5432312_SAM24385455 LIB5432313_SAM24385456 \
LIB5432314_SAM24385457 LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
mkdir -p sicelore_results/${sample_id}
java $java_xmx_flag -jar ../software/git/sicelore/Jar/Sicelore-2.0.jar \
  IsoformMatrix \
  I=bam/${sample_id}_sicelore_GE.bam \
  OUTDIR=sicelore_results/${sample_id} \
  PREFIX=sicelore \
  ISOBAM=true \
  GENETAG=GE UMITAG=U8 CELLTAG=BC \
  REFFLAT=reference_data/ensembl_90.refFlat.txt \
  CSV=reference_data/barcodes.csv \
  DELTA=2 MAXCLIP=150 METHOD=STRICT AMBIGUOUS_ASSIGN=false
done

# Align Nanopore reads to the reference transcriptome with minimap2 for NanoCount
for sample_id in LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
../software/bin/minimap2 -t $ncpu \
  -ax map-ont -p 0 -N 10 \
  ../reference_data/transcriptome.fasta \
  fastq/${sample_id}.fastq.gz \
  | samtools view -bh > bam/${sample_id}_nanocount.bam
done

# Run NanoCount
conda activate isosceles_nanocount
mkdir -p nanocount_results
for sample_id in LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
NanoCount \
  -i bam/${sample_id}_nanocount.bam \
  -o nanocount_results/${sample_id}.tsv
done
conda deactivate

# Run Flair
conda activate isosceles_flair
for sample_id in LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
mkdir -p flair_results/${sample_id}
echo sample$'\t'condition$'\t'batch$'\t'fastq/${sample_id}.fastq.gz \
  > flair_results/${sample_id}/read_manifest.tsv
flair 1234 \
  --threads $ncpu \
  --reads fastq/${sample_id}.fastq.gz \
  --genome  ../reference_data/genome.fasta \
  --gtf ../reference_data/ensembl_90.gtf \
  --junction_bed ../reference_data/known_introns.bed \
  --reads_manifest flair_results/${sample_id}/read_manifest.tsv \
  --output flair_results/${sample_id}/flair
done
conda deactivate

# Prepare the reference annotations (refgene) file for LIQA
conda activate isosceles_liqa
mkdir -p reference_data
liqa -task refgene -format gtf \
  -ref ../reference_data/ensembl_90.gtf \
  -out reference_data/ensembl_90.refgene
conda deactivate

# Filter the BAM files for LIQA
for sample_id in LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
samtools view \
  bam/${sample_id}.bam \
  -F 2308 -q 50 -O BAM \
  -o bam/${sample_id}_liqa.bam
samtools index bam/${sample_id}_liqa.bam
done

# Run LIQA
conda activate isosceles_liqa
mkdir -p liqa_results
for sample_id in LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
liqa -task quantify \
  -max_distance 10 -f_weight 1 \
  -refgene reference_data/ensembl_90.refgene \
  -bam bam/${sample_id}_liqa.bam \
  -out liqa_results/${sample_id}.tsv
done
conda deactivate

# Run ESPRESSO_S
conda activate isosceles_espresso_deps
for sample_id in LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
mkdir -p espresso_results/${sample_id}
echo bam/${sample_id}.bam$'\t'sample \
  > espresso_results/${sample_id}/samples.tsv
perl ../software/git/espresso/src/ESPRESSO_S.pl \
  -L espresso_results/${sample_id}/samples.tsv \
  -F ../reference_data/genome.fasta \
  -A ../reference_data/ensembl_90.gtf \
  -M MT \
  -T $ncpu \
  -O espresso_results/${sample_id}
done
conda deactivate

# Run ESPRESSO_C
conda activate isosceles_espresso_deps
for sample_id in LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
perl ../software/git/espresso/src/ESPRESSO_C.pl \
  -I espresso_results/${sample_id} \
  -F ../reference_data/genome.fasta \
  -T $ncpu \
  -X 0
done
conda deactivate

# Run ESPRESSO_Q (transcript quantification)
conda activate isosceles_espresso_deps
for sample_id in LIB5432315_SAM24385458 LIB5432316_SAM24385459 \
LIB5427896_SAM24376275 LIB5427897_SAM24376276; do
perl ../software/git/espresso/src/ESPRESSO_Q.pl \
  -L espresso_results/${sample_id}/samples.tsv.updated \
  -A ../reference_data/ensembl_90.gtf \
  -T $ncpu
done
conda deactivate
