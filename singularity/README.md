# Preparing the Singularity images

## Content

  * **run.sh** - the script for downloading the Singularity SIF images
  * **build_images.sh** - a script containing the commands used for building the Singularity images
  * **reports_r.def** - a Singularity definition file for the reports_r.sif image
  * **reports_python.def** - a Singularity definition file for the reports_python.sif image
  * **isosceles_jupyter.yml** - a conda environment YAML file used by the reports_python.def file
  * **isosceles.def** - a Singularity definition file for the isosceles.sif image
  * **bambu.def** - a Singularity definition file for the bambu.sif image

## Singularity image list

  * **reports_r.sif** - an environment for building the RMarkdown benchmark reports
  * **reports_python.sif** - an environment for building the Jupyter benchmark reports
  * **isosceles.sif** - an environment for running Isosceles
  * **bambu.sif** - an environment for running bambu

## Required programs and libraries

  * Singularity (>= 3.5.1)
