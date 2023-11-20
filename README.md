# Isosceles paper

This repository contains scripts and data analysis for the Isosceles paper, enabling the replication of data, analyses, and plots featured in the paper. The repository is divided into two main sections:

## Section #1. Analysis of Quantifications

In this section, you can access singularity images and benchmark program outputs required to run our scripts in the `reports` folder, reproducing the figures and plots presented in the paper. Each directory includes a script to download the results from running all programs on the corresponding benchmark. You can execute these scripts using the provided commands below. For more details on the exact commands used to generate the data, please refer to Section #2 or check the brief summary in the [Benchmark_commands.md](Benchmark_commands.md) document.

Download and prepare the Singularity images:
```bash
cd singularity
bash run.sh
cd ..
```

Prepare reference genome and annotations data:
```bash
cd reference_data
bash run.sh
cd ..
```

Download the data needed to create the reports:
```bash
cd simulated_bulk_analysis
bash download_report_data.sh
cd ..
cd simulated_sc_analysis
bash download_report_data.sh
cd ..
cd nanopore_bulk_analysis
bash download_report_data.sh
cd ..
cd nanopore_sc_analysis
bash download_report_data.sh
cd ..
```

Then you can run the desired reports to generate figures - this can be done manually in RStudio and Jupyter, or using the following commands (the latter requires Singularity to be installed):
```bash
cd reports
bash run.sh
cd ..
```

## Section #2. Mouse E18 Brain scRNA-Seq (Lebrigand et al.) Data Analysis

The analysis of the Lebrigand et al. data from the paper can be replicated using the code in the `mouse_E18_analysis` directory. You can download the input data for the analysis (including scRNA-Seq BAM files prepared using the Sicelore workflow described in Part #1), or start from any downstream step, as RDS files containing data necessary for further analysis are provided for every script.

Download and prepare the Singularity images:
```bash
cd singularity
bash run.sh
cd ..
```

Prepare input data for the analysis:
```bash
cd mouse_E18_analysis
bash download_data.sh
cd ..
```

The `mouse_E18_analysis` directory contains the following scripts for generating the figures and tables:
```bash
00_run_isosceles.R
01_scrnaseq_analysis.R
02_run_dexseq_intra.R
03_run_dexseq_inter.R
04_psi_heatmap_intra.R
05_psi_heatmap_inter.R
06_other_plots.R
07_isoswitch_analysis.R
```

You can run these scripts manually in RStudio, or using the following commands (the latter requires Singularity to be installed):
```bash
cd mouse_E18_analysis
bash run.sh
cd ..
```


## Supplementary Section. (Optional) Re-generating Quantifications from Raw Reads

We also provide the scripts to re-run all the programs to re-generate quantifications used in Part #1 `reports` (provided you download the input data).  This is more involved, and requires significant computation and memory resources.  While Isosceles is fairly efficient and can be run on most systems, not all software benchmarked here is, and so most steps were run on a cluster with 20 CPUs and 200 GB of RAM.  While we provide environments for each software suite, we note that although this runs cleanly on our system, some programs may take some troubleshooting on other systems to get working.

First prepare the software and environments using:
```bash
cd software
bash run.sh
```

Prepare the Singularity images:
```bash
cd singularity
bash run.sh
cd ..
```

Prepare reference genome and annotations data using:
```bash
cd reference_data
bash run.sh
cd ..
```

Prepare input data for the analysis using (Note: this script also contains 
information about manually downloading data from EGA):
```bash
cd input_data
bash run.sh
cd ..
```

Then in each of the following directories you can follow the `run.sh` steps:
```bash
cd simulated_bulk_analysis
bash run.sh
cd ..
cd simulated_sc_analysis
bash run.sh
cd ..
cd illumina_sc_analysis
bash run.sh
cd ..
cd nanopore_bulk_analysis
bash run.sh
cd ..
cd nanopore_sc_analysis
bash run.sh
cd ..
```

Finally to prepare the data for the analysis scripts in `reports`, as used in Part #1, you can run the `prepare_report_data.sh` (or `prepare_report_data.R`, depending on the directory) scripts:
```bash
cd simulated_bulk_analysis
singularity exec ../singularity/isosceles.sif Rscript prepare_report_data.R
cd ..
cd simulated_sc_analysis
bash prepare_report_data.sh
cd ..
cd nanopore_bulk_analysis
bash prepare_report_data.sh
cd ..
cd nanopore_sc_analysis
bash prepare_report_data.sh
cd ..
```
