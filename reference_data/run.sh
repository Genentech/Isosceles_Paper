#! /usr/bin/env bash

# Download reference genome and annotation data
wget https://zenodo.org/record/10910103/files/reference_data.tgz
tar xfz reference_data.tgz
rm -f reference_data.tgz
