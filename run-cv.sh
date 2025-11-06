#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=new-cv
#SBATCH --constraint="cpuonly"
#SBATCH --output=out/lacv_%A.out
#SBATCH --error=err/lacv_%A.err
#SBATCH --ntasks=80
#SBATCH --nodes=1
#SBATCH --mem=60GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ldacunha@ucsc.edu

set -euo pipefail

module load gsl/2.8 udunits/2.2.28-gcc14.2
module load freetype/2.12.1
module load gdal/3.9.2 r/4.4.1

export K=10
export SAMPLES=1000
export WARMUP=1000
export CHAINS=4
export CORES=1
export DEBUG=0

Rscript code/la-10-par.R
