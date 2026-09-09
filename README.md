# VIPER (VIral CaPture sEquencing Runner)
Snakemake Workflow for VirCapSeq <br>
<img width="357" height="827" alt="pipeline_graph" src="https://github.com/user-attachments/assets/ef01a88a-1317-41f2-94ae-1eb9a9cef020" />


# Installation
To install the pipeline a conda environment will be created by using the environment.yaml file, which will install all software necessary for this pipeline, make sure `environment.yaml` `Snakefile` `snakemake_script.sh` are in the same directory <br>

## Install and activate the conda environment
`conda env create --name snakemake_vircapseq --file environment.yaml` <br>
`conda activate snakemake_vircapseq` <br>

## Install the pandas module for the pipeline
`pip install pandas` <br>

# Configuration
To configure certain parameters and file paths of your choice, the `config.yaml` file will need to be adjusted by replacing some of the values. The source code (`Snakefile`) may be edited for the number of threads for each rule (will be changed in a later update) <br>

# Submitting to a job scheduler (HPC)
To submit the pipeline to your HPC, a shell script `snakemake_script.sh` has been provided for you, only adjustments need to be made are specifying the number of `--cores` available on your HPC and `SBATCH` requests <br>

# Running the pipeline
Before running the pipeline do: <br>
`snakemake -np` <br>

This does a dry run of the pipeline to make sure all files are connected and the number of files generated at each table. If no error message shows up, then you are ready to submit <br>
`sbatch snakemake_script.sh`
