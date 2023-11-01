#! /usr/bin/env bash

# Build the Singularity SIF images
sudo singularity build isosceles.sif isosceles.def
sudo singularity build bambu.sif bambu.def
sudo singularity build reports_python.sif reports_python.def
