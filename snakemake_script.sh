#!/bin/bash
#SBATCH --job-name=snakemake_bats
#SBATCH -A zedru
#SBATCH -p ceres
#SBATCH -N 1
#SBATCH -n 96
#SBATCH --mem=2000G
#SBATCH -t 300:00:00
#SBATCH --mail-user=bryan.zavalamartinez@usda.gov
#SBATCH --mail-type=BEGIN,END,FINAL
#SBATCH -o "stdout/stdout.snake.%j"
#SBATCH -e "stderr/stderr.snake.%j"

module load miniconda
source activate snakemake_rodent_data

snakemake --cores 2000 --latency-wait 60 --rerun-incomplete --nolock 
