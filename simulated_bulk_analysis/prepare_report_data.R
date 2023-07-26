#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(Isosceles)

# Global parameters
data_dir <- "report_data"
dir.create(data_dir, recursive = TRUE)
sample_id <- "truncated_bulk_rnaseq"
perc_downs <- c(10, 20, 30)

# Prepare reference data
dir.create(file.path(data_dir, "reference"), recursive = TRUE)
file.copy(file.path("..", "reference_data", "benchmark_transcript_annotations.gtf"),
          file.path(data_dir, "reference", "benchmark_full_annotations.gtf"),
          overwrite = TRUE)
file.copy(file.path("..", "reference_data", "benchmark_downsampled_10.gtf"),
          file.path(data_dir, "reference", "benchmark_downsampled_10.gtf"),
          overwrite = TRUE)
file.copy(file.path("..", "reference_data", "benchmark_downsampled_20.gtf"),
          file.path(data_dir, "reference", "benchmark_downsampled_20.gtf"),
          overwrite = TRUE)
file.copy(file.path("..", "reference_data", "benchmark_downsampled_30.gtf"),
          file.path(data_dir, "reference", "benchmark_downsampled_30.gtf"),
          overwrite = TRUE)
file.copy(file.path("..", "reference_data", "benchmark_transcript_expression.tab"),
          file.path(data_dir, "reference", "benchmark_transcript_expression.tsv"),
          overwrite = TRUE)

# Create output data directories
dir.create(file.path(data_dir, glue("{sample_id}_quant")), recursive = TRUE)
dir.create(file.path(data_dir, glue("{sample_id}_denovo")), recursive = TRUE)
dir.create(file.path(data_dir, glue("{sample_id}_stringtie")), recursive = TRUE)


# Prepare Isosceles data (transcript quantification)
se <- readRDS(file.path("isosceles_results", glue("{sample_id}_quant_se_transcript.rds")))
tpm_df <- data.frame(
    transcript_id = rownames(se),
    tpm = assay(se, "tpm")[, 1]
)
write.table(
    tpm_df,
    file.path(data_dir, glue("{sample_id}_quant"), "isosceles.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)
export_gtf(
    se,
    file.path(data_dir, glue("{sample_id}_quant"), "isosceles.gtf")
)


# Prepare Isosceles data (de novo detection)
for (perc_down in perc_downs) {
    se <- readRDS(file.path("isosceles_results", glue("{sample_id}_denovo_{perc_down}_se_transcript.rds")))
    tx_selector <- rowData(se)$splicing_support_level == "AP" |
        assay(se, "relative_expression")[, 1] >= 0
    se <- se[tx_selector,]
    tpm_df <- data.frame(
        transcript_id = rownames(se),
        tpm = assay(se, "tpm")[, 1]
    )
    write.table(
        tpm_df,
        file.path(data_dir, glue("{sample_id}_denovo"), glue("isosceles_down{perc_down}.tsv")),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    write.table(
        tpm_df,
        file.path(data_dir, glue("{sample_id}_stringtie"), glue("isosceles_down{perc_down}.tsv")),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    export_gtf(
        se,
        file.path(data_dir, glue("{sample_id}_denovo"), glue("isosceles_down{perc_down}.gtf"))
    )
    export_gtf(
        se,
        file.path(data_dir, glue("{sample_id}_stringtie"), glue("isosceles_down{perc_down}.gtf"))
    )
}

# Prepare Isosceles data (de novo detection using StringTie)
for (perc_down in perc_downs) {
    se <- readRDS(file.path("isosceles_results", glue("{sample_id}_stringtie_{perc_down}_se_transcript.rds")))
    tx_selector <- rowData(se)$splicing_support_level == "AP" |
        assay(se, "relative_expression")[, 1] >= 0
    se <- se[tx_selector,]
    tpm_df <- data.frame(
        transcript_id = rownames(se),
        tpm = assay(se, "tpm")[, 1]
    )
    write.table(
        tpm_df,
        file.path(data_dir, glue("{sample_id}_stringtie"), glue("isosceles_stringtie_down{perc_down}.tsv")),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    export_gtf(
        se,
        file.path(data_dir, glue("{sample_id}_stringtie"), glue("isosceles_stringtie_down{perc_down}.gtf"))
    )
}

# Prepare Isosceles data (de novo detection using IsoQuant)
for (perc_down in perc_downs) {
    se <- readRDS(file.path("isosceles_results", glue("{sample_id}_isoquant_{perc_down}_se_transcript.rds")))
    tx_selector <- rowData(se)$splicing_support_level == "AP" |
        assay(se, "relative_expression")[, 1] >= 0
    se <- se[tx_selector,]
    tpm_df <- data.frame(
        transcript_id = rownames(se),
        tpm = assay(se, "tpm")[, 1]
    )
    write.table(
        tpm_df,
        file.path(data_dir, glue("{sample_id}_stringtie"), glue("isosceles_IsoQuant_down{perc_down}.tsv")),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    export_gtf(
        se,
        file.path(data_dir, glue("{sample_id}_stringtie"), glue("isosceles_IsoQuant_down{perc_down}.gtf"))
    )
}

# Prepare bambu data (transcript quantification)
tpm_df <- read.delim(file.path("bambu_results", glue("{sample_id}_quant_counts_transcript.txt")))
colnames(tpm_df) <- c("TXNAME", "GENEID", "Counts")
tpm_df <- tpm_df %>%
    transmute(
        transcript_id = TXNAME,
        tpm = Counts / sum(Counts) * 1e6
    )
write.table(
    tpm_df,
    file.path(data_dir, glue("{sample_id}_quant"), "bambu.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)
file.copy(file.path("bambu_results", glue("{sample_id}_quant_extended_annotations.gtf")),
          file.path(data_dir, glue("{sample_id}_quant"), "bambu.gtf"),
          overwrite = TRUE)

# Prepare bambu data (de novo detection)
for (perc_down in perc_downs) {
    tpm_df <- read.delim(file.path("bambu_results", glue("{sample_id}_denovo_{perc_down}_counts_transcript.txt")))
    colnames(tpm_df) <- c("TXNAME", "GENEID", "Counts")
    tpm_df <- tpm_df %>%
        transmute(
            transcript_id = TXNAME,
            tpm = Counts / sum(Counts) * 1e6
        )
    write.table(
        tpm_df,
        file.path(data_dir, glue("{sample_id}_denovo"), glue("bambu_down{perc_down}.tsv")),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    file.copy(file.path("bambu_results", glue("{sample_id}_denovo_{perc_down}_extended_annotations.gtf")),
              file.path(data_dir, glue("{sample_id}_denovo"), glue("bambu_down{perc_down}.gtf")),
              overwrite = TRUE)
    write.table(
        tpm_df,
        file.path(data_dir, glue("{sample_id}_stringtie"), glue("bambu_down{perc_down}.tsv")),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    file.copy(file.path("bambu_results", glue("{sample_id}_denovo_{perc_down}_extended_annotations.gtf")),
              file.path(data_dir, glue("{sample_id}_stringtie"), glue("bambu_down{perc_down}.gtf")),
              overwrite = TRUE)
}

# Prepare bambu data (de novo detection using StringTie)
for (perc_down in perc_downs) {
    tpm_df <- read.delim(file.path("bambu_results", glue("{sample_id}_stringtie_{perc_down}_counts_transcript.txt")))
    colnames(tpm_df) <- c("TXNAME", "GENEID", "Counts")
    tpm_df <- tpm_df %>%
        transmute(
            transcript_id = TXNAME,
            tpm = Counts / sum(Counts) * 1e6
        )
    write.table(
        tpm_df,
        file.path(data_dir, glue("{sample_id}_stringtie"), glue("bambu_stringtie_down{perc_down}.tsv")),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    file.copy(file.path("bambu_results", glue("{sample_id}_stringtie_{perc_down}_extended_annotations.gtf")),
              file.path(data_dir, glue("{sample_id}_stringtie"), glue("bambu_stringtie_down{perc_down}.gtf")),
              overwrite = TRUE)
}

# Prepare Flair data (transcript quantification)
tpm_df <- read.delim(file.path("flair_results",
                               glue("{sample_id}_quant"),
                               "flair.counts_matrix.tsv.counts.tsv"))
colnames(tpm_df) <- c("TXNAME", "Counts")
tpm_df <- tpm_df %>%
    transmute(
        transcript_id = TXNAME,
        tpm = Counts / sum(Counts) * 1e6
    )
write.table(
    tpm_df,
    file.path(data_dir, glue("{sample_id}_quant"), "flair.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)
file.copy(file.path("flair_results", glue("{sample_id}_quant"),
                    "flair.isoforms.gtf"),
          file.path(data_dir, glue("{sample_id}_quant"), "flair.gtf"),
          overwrite = TRUE)

# Prepare Flair data (de novo detection)
for (perc_down in perc_downs) {
    tpm_df <- read.delim(file.path("flair_results",
                                   glue("{sample_id}_denovo_{perc_down}"),
                                   "flair.counts_matrix.tsv.counts.tsv"))
    colnames(tpm_df) <- c("TXNAME", "Counts")
    tpm_df <- tpm_df %>%
        transmute(
            transcript_id = TXNAME,
            tpm = Counts / sum(Counts) * 1e6
        )
    write.table(
        tpm_df,
        file.path(data_dir, glue("{sample_id}_denovo"), glue("flair_down{perc_down}.tsv")),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    write.table(
        tpm_df,
        file.path(data_dir, glue("{sample_id}_stringtie"), glue("flair_down{perc_down}.tsv")),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    file.copy(file.path("flair_results", glue("{sample_id}_denovo_{perc_down}"),
                        "flair.isoforms.gtf"),
              file.path(data_dir, glue("{sample_id}_denovo"), glue("flair_down{perc_down}.gtf")),
              overwrite = TRUE)
    file.copy(file.path("flair_results", glue("{sample_id}_denovo_{perc_down}"),
                        "flair.isoforms.gtf"),
              file.path(data_dir, glue("{sample_id}_stringtie"), glue("flair_down{perc_down}.gtf")),
              overwrite = TRUE)
}

# Prepare Flair data (de novo detection using StringTie)
for (perc_down in perc_downs) {
    tpm_df <- read.delim(file.path("flair_results",
                                   glue("{sample_id}_stringtie_{perc_down}"),
                                   "flair.counts_matrix.tsv.counts.tsv"))
    colnames(tpm_df) <- c("TXNAME", "Counts")
    tpm_df <- tpm_df %>%
        transmute(
            transcript_id = TXNAME,
            tpm = Counts / sum(Counts) * 1e6
        )
    write.table(
        tpm_df,
        file.path(data_dir, glue("{sample_id}_stringtie"), glue("flair_stringtie_down{perc_down}.tsv")),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    file.copy(file.path("flair_results", glue("{sample_id}_stringtie_{perc_down}"),
                        "flair.isoforms.gtf"),
              file.path(data_dir, glue("{sample_id}_stringtie"), glue("flair_stringtie_down{perc_down}.gtf")),
              overwrite = TRUE)
}

# Prepare IsoQuant data (transcript quantification)
file.copy(file.path("isoquant_results", glue("{sample_id}_quant"),
                    glue("00_{sample_id}"),
                    glue("00_{sample_id}.transcript_tpm.tsv")),
          file.path(data_dir, glue("{sample_id}_quant"), "isoquant.tsv"),
          overwrite = TRUE)
file.copy(file.path("..", "reference_data", "benchmark_transcript_annotations.gtf"),
          file.path(data_dir, glue("{sample_id}_quant"), "isoquant.gtf"),
          overwrite = TRUE)

# Prepare IsoQuant data (de novo detection)
for (perc_down in perc_downs) {
    file.copy(file.path("isoquant_results", glue("{sample_id}_quant_denovo_{perc_down}"),
                        glue("00_{sample_id}"),
                        glue("00_{sample_id}.transcript_tpm.tsv")),
              file.path(data_dir, glue("{sample_id}_denovo"), glue("isoquant_down{perc_down}.tsv")),
              overwrite = TRUE)
    file.copy(file.path("isoquant_results", glue("{sample_id}_quant_denovo_{perc_down}"),
                        glue("00_{sample_id}"),
                        glue("00_{sample_id}.transcript_tpm.tsv")),
              file.path(data_dir, glue("{sample_id}_stringtie"), glue("isoquant_down{perc_down}.tsv")),
              overwrite = TRUE)
    file.copy(file.path("isoquant_results", glue("{sample_id}_denovo_{perc_down}"),
                        "extended_annotations.gtf"),
              file.path(data_dir, glue("{sample_id}_denovo"), glue("isoquant_down{perc_down}.gtf")),
              overwrite = TRUE)
    file.copy(file.path("isoquant_results", glue("{sample_id}_denovo_{perc_down}"),
                        "extended_annotations.gtf"),
              file.path(data_dir, glue("{sample_id}_stringtie"), glue("isoquant_down{perc_down}.gtf")),
              overwrite = TRUE)
}

# Prepare IsoQuant data (de novo detection using StringTie)
for (sample_id in sample_ids) {
    for (perc_down in perc_downs) {
        file.copy(file.path("isoquant_results", glue("{sample_id}_quant_stringtie_{perc_down}"),
                            glue("00_{sample_id}"),
                            glue("00_{sample_id}.transcript_tpm.tsv")),
                  file.path(data_dir, glue("{sample_id}_stringtie"), glue("isoquant_stringtie_down{perc_down}.tsv")),
                  overwrite = TRUE)
        file.copy(file.path("isoquant_results", glue("{sample_id}_stringtie_{perc_down}"),
                            "extended_annotations.gtf"),
                  file.path(data_dir, glue("{sample_id}_stringtie"), glue("isoquant_stringtie_down{perc_down}.gtf")),
                  overwrite = TRUE)
    }
}

# Prepare NanoCount data (transcript quantification)
file.copy(file.path("nanocount_results", glue("{sample_id}_quant.tsv")),
          file.path(data_dir, glue("{sample_id}_quant"), "nanocount.tsv"),
          overwrite = TRUE)
file.copy(file.path("..", "reference_data", "benchmark_transcript_annotations.gtf"),
          file.path(data_dir, glue("{sample_id}_quant"), "nanocount.gtf"),
          overwrite = TRUE)

# Prepare NanoCount data (de novo detection using StringTie)
for (perc_down in perc_downs) {
    file.copy(file.path("nanocount_results", glue("{sample_id}_stringtie_{perc_down}.tsv")),
              file.path(data_dir, glue("{sample_id}_stringtie"), glue("nanocount_stringtie_down{perc_down}.tsv")),
              overwrite = TRUE)
    file.copy(file.path("stringtie_results", glue("{sample_id}_denovo_{perc_down}.gtf")),
              file.path(data_dir, glue("{sample_id}_stringtie"), glue("nanocount_stringtie_down{perc_down}.gtf")),
              overwrite = TRUE)
}

# Prepare LIQA data (transcript quantification)
tpm_df <- read.delim(file.path("liqa_results", glue("{sample_id}_quant.tsv")))
tpm_df <- tpm_df %>%
    transmute(
        transcript_id = IsoformName,
        tpm = ReadPerGene_corrected / sum(ReadPerGene_corrected) * 1e6
    )
write.table(
    tpm_df,
    file.path(data_dir, glue("{sample_id}_quant"), "liqa.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)
file.copy(file.path("..", "reference_data", "benchmark_transcript_annotations.gtf"),
          file.path(data_dir, glue("{sample_id}_quant"), "liqa.gtf"),
          overwrite = TRUE)

# Prepare LIQA data (de novo detection using StringTie)
for (perc_down in perc_downs) {
    tpm_df <- read.delim(file.path("liqa_results",
                                   glue("{sample_id}_stringtie_{perc_down}.tsv")))
    tpm_df <- tpm_df %>%
        transmute(
            transcript_id = IsoformName,
            tpm = ReadPerGene_corrected / sum(ReadPerGene_corrected) * 1e6
        )
    write.table(
        tpm_df,
        file.path(data_dir, glue("{sample_id}_stringtie"),
                  glue("liqa_stringtie_down{perc_down}.tsv")),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    file.copy(file.path("stringtie_results", glue("{sample_id}_denovo_{perc_down}.gtf")),
              file.path(data_dir, glue("{sample_id}_stringtie"), glue("liqa_stringtie_down{perc_down}.gtf")),
              overwrite = TRUE)
}

# Prepare ESPRESSO data (transcript quantification)
tpm_df <- read.delim(file.path("espresso_results",
                               glue("{sample_id}_quant"),
                               "samples_N2_R0_abundance.esp"))
tpm_df <- tpm_df %>%
    transmute(
        transcript_id = transcript_ID,
        tpm = sample / sum(sample) * 1e6
    )
write.table(
    tpm_df,
    file.path(data_dir, glue("{sample_id}_quant"), "espresso.tsv"),
    row.names = FALSE, quote = FALSE, sep = "\t"
)
file.copy(file.path("espresso_results",
                    glue("{sample_id}_quant"),
                    "samples_N2_R0_updated.gtf"),
          file.path(data_dir, glue("{sample_id}_quant"), "espresso.gtf"),
          overwrite = TRUE)

# Prepare ESPRESSO data (de novo detection)
for (perc_down in perc_downs) {
    tpm_df <- read.delim(file.path("espresso_results",
                                   glue("{sample_id}_denovo_{perc_down}"),
                                   "samples_N2_R0_abundance.esp"))
    tpm_df <- tpm_df %>%
        transmute(
            transcript_id = transcript_ID,
            tpm = sample / sum(sample) * 1e6
        )
    write.table(
        tpm_df,
        file.path(data_dir, glue("{sample_id}_denovo"), glue("espresso_down{perc_down}.tsv")),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    write.table(
        tpm_df,
        file.path(data_dir, glue("{sample_id}_stringtie"), glue("espresso_down{perc_down}.tsv")),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    file.copy(file.path("espresso_results", glue("{sample_id}_denovo_{perc_down}"),
                        "samples_N2_R0_updated.gtf"),
              file.path(data_dir, glue("{sample_id}_denovo"), glue("espresso_down{perc_down}.gtf")),
              overwrite = TRUE)
    file.copy(file.path("espresso_results", glue("{sample_id}_denovo_{perc_down}"),
                        "samples_N2_R0_updated.gtf"),
              file.path(data_dir, glue("{sample_id}_stringtie"), glue("espresso_down{perc_down}.gtf")),
              overwrite = TRUE)
}

# Prepare ESPRESSO data (de novo detection using StringTie)
for (perc_down in perc_downs) {
    tpm_df <- read.delim(file.path("espresso_results",
                                   glue("{sample_id}_stringtie_{perc_down}"),
                                   "samples_N2_R0_abundance.esp"))
    tpm_df <- tpm_df %>%
        transmute(
            transcript_id = transcript_ID,
            tpm = sample / sum(sample) * 1e6
        )
    write.table(
        tpm_df,
        file.path(data_dir, glue("{sample_id}_stringtie"), glue("espresso_stringtie_down{perc_down}.tsv")),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    file.copy(file.path("espresso_results", glue("{sample_id}_stringtie_{perc_down}"),
                        "samples_N2_R0_updated.gtf"),
              file.path(data_dir, glue("{sample_id}_stringtie"), glue("espresso_stringtie_down{perc_down}.gtf")),
              overwrite = TRUE)
}

# Print session information
sessionInfo()
