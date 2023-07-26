#! /usr/bin/env bash

# Set the number of CPUs/threads for the analysis
export ncpu=1

# Set the Java max memory flag
export java_xmx_flag="-Xmx1g"

# Set the temporary directory path
export temp_dir="/tmp/"

# Set the sample ID
export sample_id="LIB5445493_SAM24404003"

# Prepare the FASTQ file
mkdir -p fastq
ln -s ../input_data/fastq_ont/${sample_id}.fastq.gz \
  fastq/${sample_id}.fastq.gz

# Process the reference annotations
mkdir -p reference_data
cat ../reference_data/ensembl_90.gtf | awk '$3 == "exon"' > \
  reference_data/ensembl_90_exon.gtf

# Process the FASTQ file with match_cell_barcode for FLAMES
conda activate isosceles_flames_deps
mkdir -p flames_match_bc
../software/git/FLAMES/src/bin/match_cell_barcode \
  fastq \
  flames_match_bc/${sample_id}_bc_stats.txt \
  flames_match_bc/${sample_id}.fastq.gz \
  ../illumina_sc_analysis/data/barcodes.csv \
  1 12
conda deactivate

# Run FLAMES
conda activate isosceles_flames_deps
mkdir -p flames_results
../software/git/FLAMES/python/sc_long_pipeline.py \
  -a reference_data/ensembl_90_exon.gtf \
  -i flames_match_bc/${sample_id}.fastq.gz \
  --outdir flames_results \
  --genomefa ../reference_data/genome.fasta \
  --minimap2_dir $(which minimap2 | xargs dirname)
conda deactivate

# Prepare the reference annotations (refFlat) file for Sicelore
mkdir -p reference_data
gtfToGenePred -genePredExt -ignoreGroupsWithoutExons \
  ../reference_data/ensembl_90.gtf \
  /dev/stdout | \
  awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > \
  reference_data/ensembl_90.refFlat.txt

# Sicelore: Nanopore poly(A) scan, stranding of reads
mkdir -p sicelore
java $java_xmx_flag -jar ../software/git/sicelore/Jar/NanoporeReadScanner-0.5.jar \
  -i fastq/${sample_id}.fastq.gz \
  -o sicelore

# Sicelore: splitting the FASTQ file into chunks
mkdir -p sicelore/fastq_chunk
../software/bin/fastp --thread $ncpu -Q -A \
  --split_prefix_digits=1  --split=8 \
  -i sicelore/passed/${sample_id}.fastqFWD.gz \
  --out1=sicelore/fastq_chunk/chunk.fastq.gz

# Sicelore: mapping of Nanopore reads to the reference genome with minimap2
mkdir -p sicelore/bam_chunk
for i in `seq 8`; do
../software/bin/minimap2 -t $ncpu \
  -ax splice -uf --MD --sam-hit-only \
  --junc-bed ../reference_data/known_introns.bed --junc-bonus 15 \
  ../reference_data/genome.fasta \
  sicelore/fastq_chunk/$i.chunk.fastq.gz \
  | samtools sort -o sicelore/bam_chunk/$i.chunk.bam
samtools index sicelore/bam_chunk/$i.chunk.bam
done

# Sicelore: tag Nanopore SAM records with gene names
for i in `seq 8`; do
java $java_xmx_flag -jar ../software/git/sicelore/Jar/Sicelore-2.0.jar \
  AddGeneNameTag \
  I=sicelore/bam_chunk/$i.chunk.bam \
  O=sicelore/bam_chunk/$i.chunk.GE.bam \
  REFFLAT=reference_data/ensembl_90.refFlat.txt \
  ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true \
  VALIDATION_STRINGENCY=SILENT
samtools index sicelore/bam_chunk/$i.chunk.GE.bam
done

# Sicelore: add read sequence and QV values to Nanopore SAM records
for i in `seq 8`; do
mkdir -p sicelore/fastq_chunk/$i
zcat sicelore/fastq_chunk/$i.chunk.fastq.gz > \
  sicelore/fastq_chunk/$i/$i.chunk.fastq
java $java_xmx_flag -jar ../software/git/sicelore/Jar/Sicelore-2.0.jar \
  AddBamReadSequenceTag \
  I=sicelore/bam_chunk/$i.chunk.GE.bam \
  O=sicelore/bam_chunk/$i.chunk.GE_US.bam \
  FASTQDIR=sicelore/fastq_chunk/$i
rm -rf sicelore/fastq_chunk/$i
samtools index sicelore/bam_chunk/$i.chunk.GE_US.bam
done

# Sicelore: barcode and UMI assignment to Nanopore SAM records
for i in `seq 8`; do
java $java_xmx_flag -jar ../software/git/sicelore/Jar/NanoporeBC_UMI_finder-1.0.jar \
  -i sicelore/bam_chunk/$i.chunk.GE_US.bam \
  -o sicelore/bam_chunk/$i.chunk.GE_US_10xAttr.bam \
  -k ../illumina_sc_analysis/data/illumina_parsed.obj \
  --ncpu $ncpu \
  --maxUMIfalseMatchPercent 1 --maxBCfalseMatchPercent 5 \
  --logFile sicelore/bam_chunk/NanoporeBC_UMI_finder.$i.log
samtools index sicelore/bam_chunk/$i.chunk.GE_US_10xAttr.bam
samtools index sicelore/bam_chunk/$i.chunk.GE_US_10xAttr_umifound_.bam
done

# Sicelore: merge chunk bam files
samtools merge \
  sicelore/minimap2_GE_US_10xAttr_umifound_.bam \
  sicelore/bam_chunk/1.chunk.GE_US_10xAttr_umifound_.bam \
  sicelore/bam_chunk/2.chunk.GE_US_10xAttr_umifound_.bam \
  sicelore/bam_chunk/3.chunk.GE_US_10xAttr_umifound_.bam \
  sicelore/bam_chunk/4.chunk.GE_US_10xAttr_umifound_.bam \
  sicelore/bam_chunk/5.chunk.GE_US_10xAttr_umifound_.bam \
  sicelore/bam_chunk/6.chunk.GE_US_10xAttr_umifound_.bam \
  sicelore/bam_chunk/7.chunk.GE_US_10xAttr_umifound_.bam \
  sicelore/bam_chunk/8.chunk.GE_US_10xAttr_umifound_.bam
samtools index sicelore/minimap2_GE_US_10xAttr_umifound_.bam

# Sicelore: split bam by chromosomes
mkdir -p sicelore/bam_split
chrom_names=$(grep ">" ../reference_data/genome.fasta | cut -d " " -f1 | sed 's/>//' -)
for chrom_name in $chrom_names; do
samtools view \
  -Sb sicelore/minimap2_GE_US_10xAttr_umifound_.bam \
  ${chrom_name} \
  -o sicelore/bam_split/minimap2_GE_US_10xAttr_umifound_${chrom_name}.bam
samtools index sicelore/bam_split/minimap2_GE_US_10xAttr_umifound_${chrom_name}.bam
done

# Sicelore: generate consensus sequences
mkdir -p sicelore/fastq_split
chrom_names=$(grep ">" ../reference_data/genome.fasta | cut -d " " -f1 | sed 's/>//' -)
for chrom_name in $chrom_names; do
java $java_xmx_flag -jar ../software/git/sicelore/Jar/Sicelore-2.0.jar \
  ComputeConsensus \
  I=sicelore/bam_split/minimap2_GE_US_10xAttr_umifound_${chrom_name}.bam \
  O=sicelore/fastq_split/molecules.${chrom_name}.fastq \
  T=$ncpu TMPDIR=$temp_dir
done

# Sicelore: deduplicate molecules from genes having multiple copies in the genome
cat sicelore/fastq_split/molecules.*.fastq > \
  sicelore/molecules.fastq
java $java_xmx_flag -jar ../software/git/sicelore/Jar/Sicelore-2.0.jar \
  DeduplicateMolecule \
  I=sicelore/molecules.fastq \
  O=sicelore/molecules_deduplicated.fastq

# Sicelore: mapping of molecules consensus sequences to the reference genome with minimap2
../software/bin/minimap2 -t $ncpu \
  -ax splice -uf --sam-hit-only --secondary=no \
  --junc-bed ../reference_data/known_introns.bed --junc-bonus 15 \
  ../reference_data/genome.fasta \
  sicelore/molecules_deduplicated.fastq \
  | samtools sort -o sicelore/molecules.bam
samtools index sicelore/molecules.bam

# Sicelore: tag molecule SAM records with gene names, cell barcodes, UMI sequence and reads number
java $java_xmx_flag -jar ../software/git/sicelore/Jar/Sicelore-2.0.jar \
  AddGeneNameTag \
  I=sicelore/molecules.bam \
  O=sicelore/molecules_GE.bam \
  REFFLAT=reference_data/ensembl_90.refFlat.txt
samtools index sicelore/molecules_GE.bam
java $java_xmx_flag -jar ../software/git/sicelore/Jar/Sicelore-2.0.jar \
  AddBamMoleculeTags \
  I=sicelore/molecules_GE.bam \
  O=sicelore/molecules_GE_tags.bam
samtools index sicelore/molecules_GE_tags.bam

# Run Sicelore
mkdir -p sicelore_results
java $java_xmx_flag -jar ../software/git/sicelore/Jar/Sicelore-2.0.jar \
  IsoformMatrix \
  I=sicelore/molecules_GE_tags.bam \
  OUTDIR=sicelore_results \
  PREFIX=sicelore \
  ISOBAM=true \
  GENETAG=GE UMITAG=U8 CELLTAG=BC \
  REFFLAT=reference_data/ensembl_90.refFlat.txt \
  CSV=../illumina_sc_analysis/data/barcodes.csv \
  DELTA=2 MAXCLIP=150 METHOD=STRICT AMBIGUOUS_ASSIGN=false

# Run IsoQuant
conda activate isosceles_isoquant
mkdir -p isoquant_results
isoquant.py \
  --threads $ncpu \
  --no_model_construction \
  --data_type nanopore \
  --reference ../reference_data/genome.fasta \
  --genedb ../reference_data/ensembl_90.gtf \
  --bam sicelore/molecules_GE_tags.bam \
  --read_group tag:BC \
  -o isoquant_results
conda deactivate

# Run Isosceles
singularity exec ../singularity/isosceles.sif Rscript run_isosceles.R
