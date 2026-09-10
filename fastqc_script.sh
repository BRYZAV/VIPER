#!/bin/bash
#SBATCH --job-name="fastqc_script"
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 10:00:00


module load miniconda
source activate metagenomics

mkdir -p Raw_reads_QC
mkdir -p Filtered_reads_QC

fastqc fastq_data/*gz -t 48 -o Raw_reads_QC/
multiqc Raw_reads_QC/ -n Raw_MultiQC_report.html -o Raw_MultiQC_report
fastqc trimmed_reads/*gz -t 48 -o Filtered_reads_QC/
multiqc Filtered_reads_QC/ -n Filtered_MultiQC_report.html -o Filtered_MultiQC_report
