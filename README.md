# Isosceles paper

Isosceles paper scripts and data analysis repository to reproduce data/analysis/plots from the paper. This repo has two parts:

## Section #1. Analysis of Quantifications

Here, you can download singularity images and the output of each benchmark program.  With these you can run our scripts in `reports` to reproduce the figures/plots in the paper. Each directory has a script to first download the results from running all programs on the that benchmark, which can be run using the commands described below. The exact commands used to generate the data are described in Section #2, but you can find a brief summary in the [Benchmark_commands.md](Benchmark_commands.md) document.

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

## Section #2. (Optional) Re-generating Quantifications from Raw Reads

We also provide the scripts to re-run all the programs to re-generate quantifications used in Part #1 `reports` (provided you download the input data).  This is quite a bit more involved, and requires significant computation and memory resources.  While Isosceles is fairly efficient and can be run on most systems, not all software benchmarked here is, and so we ran most steps on our cluster with 20 CPUs and 200 GB of RAM.  While we provide environments for each software suite, we also provide the disclaimer that although this runs cleanly on our system, some programs may take some troubleshooting on other systems to get working.

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
