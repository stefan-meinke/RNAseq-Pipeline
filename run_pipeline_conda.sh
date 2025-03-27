#!/bin/bash
# run_pipeline.sh - Run the RNAseq pipeline inside a Conda environment

# Initialize conda (adjust the path if necessary)
source "$(conda info --base)/etc/profile.d/conda.sh"

# Activate your conda environment
conda activate rnaseq_pipeline

# Use Guix shell with manifest.scm to create an isolated environment
snakemake -s pipeline_deseq_guix_shell.smk --profile slurm --cluster-config cluster.yml \
    --cluster "sbatch --mem={cluster.mem} --time={cluster.time} --cpus-per-task={cluster.cpus-per-task} {cluster.log}" \
    -j 100 > progress.log 2>&1
