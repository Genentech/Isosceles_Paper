#! /usr/bin/env bash

# Build software from source
## minimap2
git clone https://github.com/lh3/minimap2
cd minimap2
git checkout 50a26a60a6586aa8d9ac314ff0502e2b2109549c
make
cp -f minimap2 ../bin
cd ..
rm -rf minimap2
## Racon
git clone --recursive https://github.com/lbcb-sci/racon.git
cd racon
git checkout a2cfcac281d312a73912a97d6d960404f516c389
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cp bin/racon ../../bin
cd ../..
rm -rf racon
## spoa
git clone https://github.com/rvaser/spoa
cd spoa
git checkout 981ad5921c989948b37e821f6f47e8e23873cf28
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cp bin/spoa* ../../bin
cd ../..
rm -rf spoa

# Download software
## fastp
wget -O bin/fastp http://opengene.org/fastp/fastp.0.23.2
chmod a+x bin/fastp

# Create conda environments from the provided YAML files
conda env create -f conda_env/isosceles_espresso_deps.yml
conda env create -f conda_env/isosceles_flair.yml
conda env create -f conda_env/isosceles_flames_deps.yml
conda env create -f conda_env/isosceles_isoquant.yml
conda env create -f conda_env/isosceles_liqa.yml
conda env create -f conda_env/isosceles_nanocount.yml
conda env create -f conda_env/isosceles_nanosim.yml
conda env create -f conda_env/isosceles_stringtie.yml
conda env create -f conda_env/isosceles_umitools.yml

# Clone program repositories from GitHub
cd git
## ESPRESSO
git clone https://github.com/Xinglab/espresso.git
cd espresso
git checkout f7289a72be8d744acd1bcbb9d3a5ab8800131886
cd ..
## FLAMES
git clone https://github.com/LuyiTian/FLAMES.git
cd FLAMES
git checkout 774e16ae53a1430e03081970827e93f1fbaecead
cd ..
## Build the match_cell_barcode binary (FLAMES)
cd FLAMES/src
g++ -std=c++11 -lz -O2 -o bin/match_cell_barcode \
  ssw/ssw_cpp.cpp ssw/ssw.c match_cell_barcode.cpp kseq.h edit_dist.cpp
cd ../..
## Sicelore
git clone https://github.com/ucagenomix/sicelore.git
cd sicelore
git checkout 1219f60af40c30ed20f74c79b380673c992326fa
cd ..
cd ..

# Download Nextflow workflows
## wf-single-cell
nextflow pull epi2me-labs/wf-single-cell -r master
