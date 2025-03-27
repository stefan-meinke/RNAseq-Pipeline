# RNAseq Pipeline

This repository contains an RNAseq analysis pipeline built with [Snakemake](https://snakemake.readthedocs.io/en/stable/). The pipeline is designed to be modular and reproducible and integrates multiple analysis steps—from read preprocessing to differential gene expression and alternative splicing analysis and exon-level PSI calculation. It is designed to run either in a GUIX or Conda environment

## Pipeline Overview
The pipeline performs the following key steps:
- **Preprocessing**:
Uses [fastp](https://github.com/OpenGene/fastp) for adapter trimming and quality filtering.
- **Quality Control**:
Runs [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and aggregates reports with [MultiQC](https://github.com/MultiQC/MultiQC).
- **Alignment**:
Aligns reads to the reference genome using [STAR](https://pmc.ncbi.nlm.nih.gov/articles/PMC3530905/) in two-pass mode.
- **Gene-level Read Counting**:
Uses [featureCounts](https://academic.oup.com/bioinformatics/article/30/7/923/232889) (part of the subread package) to generate gene-level counts for differential expression analysis.
- **Transcript Quantification**:
Quantifies transcript abundance with [Salmon](https://combine-lab.github.io/salmon/getting_started/).
- **Alternative Splicing Analysis**:
Performs alternative splicing analysis using [rMATS-turbo](https://github.com/Xinglab/rmats-turbo) and [SPLASH2](https://github.com/refresh-bio/SPLASH/wiki).
- **PSI Calculation**:
The Python scripts used to calculate exon-level Percent Spliced In (PSI) values are adapted from [Schafer et al. 2015](https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/0471142905.hg1116s87) and originally available at https://github.com/MIAOKUI/PSI. **Note**: These scripts are distributed under the GNU General Public License v3 and are not provided as part of this repository. For integration into this pipeline see [Prerequisites](#prerequisites).

- **SPLASH2**
Runs [SPLASH2](https://github.com/refresh-bio/SPLASH) ([Chaung et al., 2023](https://www.cell.com/cell/fulltext/S0092-8674(23)01179-0), [Kokot et al., 2024](https://www.cell.com/cell/fulltext/S0092-8674(23)01179-0)) for unsupervised and reference-free k-mers of anchor and targte sequences.

## Prerequisites
1. Before running the pipeline, ensure that the **STAR and Salmon indexes** are available. You have two options:
    1. **Build the indexes** as part of the pipeline (recommended if you do not have pre-built indexes).
    2. **Reuse pre-built indexes** by setting the appropriate configuration flags in the [config.yml](./config.yml) file.
    3. **Build the indexes** outside of the pipeline using the following: 
\
        **STAR Index**:
        Generate a STAR index from your reference genome. For example run:

        ```bash
        STAR --runThreadN 16 --runMode genomeGenerate \
        --genomeDir path/to/output/index \
        --genomeFastaFiles path/to/genome/annotation.fa \
        --sjdbGTFfile path/to/gene/annotation.gtf \
        --sjdbOverhang 149 \ # general rule: read length -1
        ```

        **Salmon Index**:
        Create a Salmon index from your reference transcriptome. For example run:
        ```bash
        salmon index -t path/to/transcriptome.fa.gz -i path/to/output
        ```

**2. Reduced GFF Annotation**:
Prepare a reduced GFF annotation file (e.g. using [dexseq_prepare_annotation.py](https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html#24_preparing_the_annotation)) for PSI calculations.

**3. Snakemake Profile**:
Set up a [Snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for your computing environment (this pipeline assumes a [Slurm profile](https://github.com/snakemake-profiles)).

**4. GUIX shell**: 
if you are using the GUIX shell approach, you need to install [rMATS-turbo](https://github.com/Xinglab/rmats-turbo) and update the path to the executable in the [config.yml](./config.yml) file

**5. PSI scripts**
Make sure to [download](https://github.com/MIAOKUI/PSI) and adapt the original Python 2 scripts to Python 3 to integrate them into this pipeline. The following code snippets illustrate the required changes:

| script             | original code                                  | new code                        |
|--------------------|------------------------------------------------|---------------------------------|
| exclusion_count.py | `fout = open(outFile+".exclusion","w")`         | `fout = open(outFile,"w")`        |
| psi_calculation.py | `fout = open(basename + '.psi', 'wt')`          | `fout = open(basename, 'wt')`     |

**5. SPLASH2**
Make sure to [install SPLASH2](https://github.com/refresh-bio/SPLASH/wiki/Installation) first and update the path to the executable in the config.yml file

## Tools & Dependencies
The pipeline’s software environment is managed using Conda or GUIX. For full details see the environment.yml or the manifest.scm file, respectively. Key tools (and their versions) include:
- Python 3.10.7
- Snakemake 7.7.0
- fastp 0.23.2
- MultiQC 1.14
- STAR 2.7.8a
- samtools 1.19
- rMATS-turbo 4.3.0
- subread 2.0.3 (featureCounts)
- Salmon 1.10.1
- bedtools 2.30.0
- R-base 4.4.2 & Bioconductor DEXSeq 1.52
- HTSeq 2.0.2
- pysam 0.20.0



## Installation
**1. Clone the Repository**:
```bash
git clone https://github.com/your_username/rnaseq_pipeline.git
cd rnaseq_pipeline
```


**2.1 Set Up the Conda Environment**:
Make sure Conda is installed, then run:

```bash
conda env create -f environment.yml
conda activate rnaseq_pipeline
```

**2.2 Set up the Guix Shell**
Make sure Guix is installed, then run the following to load the Guix environment:

```bash
guix shell -m manifest.scm
```

Both 2.1 and 2.2 are not necessary if you want to run whole pipeline as is, as the run_pipeline.sh already includes loading the respective environment.

**3. Prepare Input Data**:
Configure your input FASTQ files, STAR index, Salmon index, and annotation files by editing the [config.yml](./config.yml).




## Running the Pipeline
To run the pipeline on a Slurm-based HPC cluster, use the provided [run_pipeline_guix.sh](./run_pipeline_guix.sh) or [run_pipeline_conda.sh](./run_pipeline_conda.sh) script, respectively. The script leverages your Snakemake profile (see [cluster.yml](./cluster.yml) for resource configurations) and submits jobs via Slurm. For example:

``` bash
./run_pipeline_guix.sh
```

The command executed by this script is:

```bash
snakemake -s pipeline.smk --profile slurm --cluster-config cluster.yml \
  --cluster "sbatch --mem={cluster.mem} --time={cluster.time} --cpus-per-task={cluster.cpus-per-task} {cluster.log}" \
  -j 100 > progress.log 2>&1
```

All progress and log information is directed to progress.log.

## Pipeline Structure
- pipeline.smk:
Main Snakemake workflow file that defines all pipeline rules.
- cluster.yml:
Cluster resource configuration for job submission (Slurm).
- config.yml:
User-configurable parameters including file paths and tool settings.
- environment.yml:
Conda environment specification listing all required tools and dependencies.
- run_pipeline.sh:
Shell script for submitting the pipeline (should be run from the project root on a login node).
- PSI_scripts:
Folder containing Python scripts for PSI calculation and related analyses. The scripts are from [Schafer et al. 2015](https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/0471142905.hg1116s87) (https://github.com/MIAOKUI/PSI)

## Output
The pipeline produces multiple output files:
- **Preprocessed Reads**:
Trimmed FASTQ files generated by fastp.
- **Quality Reports**:
Individual FastQC reports and a combined MultiQC report.
- **Alignment Files**:
STAR alignment outputs.
- **Gene Counts**:
Gene-level read count files from featureCounts.
- **Transcript Quantification**:
Salmon quantification outputs.
- **Alternative Splicing**:
Splicing analysis results from rMATS-turbo.
- **PSI Values**:
Exon-level PSI values the python scripts from https://github.com/MIAOKUI/PSI.
- **SPLASH2**:
Unsupervised and reference-free k-mers of anchor and target sequences, e.g. for downstream splicing analysis (see: https://github.com/refresh-bio/SPLASH)

## Troubleshooting
- Resource Issues:
If jobs fail due to memory errors, adjust resource settings in cluster.yml or your Snakemake rules.
- Input Paths:
Ensure that all file paths in config.yml are correct and accessible.
- Environment Problems:
Verify that the Conda environment is activated and that all dependencies are correctly installed.

## License
This project is licensed under the GNU General Public License v3.0. See the LICENSE file for details.
