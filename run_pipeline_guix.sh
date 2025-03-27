#!/bin/bash
# run_pipeline.sh - Run the RNAseq pipeline inside a Guix environment

# Use Guix shell with manifest.scm to create an isolated environment
guix shell -m manifest.scm -- \
    snakemake -s pipeline_guix.smk \
        --profile slurm \
        --cluster-config cluster.yml \
        --cluster "sbatch --mem={cluster.mem} --time={cluster.time} --cpus-per-task={cluster.cpus-per-task} {cluster.log}" \
        -j 100 > progress.log 2>&1
