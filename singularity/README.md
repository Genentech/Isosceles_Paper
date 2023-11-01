# Preparing the Singularity images

## Content

  * **run.sh** - the script for downloading the Singularity SIF images
  * **build_images.sh** - a script containing the commands used for building the Singularity images
  * **isosceles.def** - a Singularity definition file for the isosceles.sif image
  * **bambu.def** - a Singularity definition file for the bambu.sif image
  * **reports_python.def** - a Singularity definition file for the reports_python.sif image
  * **isosceles_jupyter.yml** - a conda environment YAML file used by the reports_python.def file

## Singularity image list

  * **isosceles.sif** - an environment for running Isosceles, building the RMarkdown reports etc.
  * **bambu.sif** - an environment for running bambu
  * **reports_python.sif** - an environment for building the Jupyter benchmark reports

## Required programs and libraries

  * Singularity (>= 3.5.1)
