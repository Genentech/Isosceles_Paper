#! /usr/bin/env bash

# Download reference genome and annotation data
wget https://zenodo.org/record/8180605/files/reference_data.tgz
tar xfz reference_data.tgz
rm -f reference_data.tgz
