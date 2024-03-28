#! /usr/bin/env bash

# Download the already prepared report data
wget https://zenodo.org/record/0000000/files/report_data_pacbio_analysis.tgz
tar xfz report_data_pacbio_analysis.tgz
rm -f report_data_pacbio_analysis.tgz
