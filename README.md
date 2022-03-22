# metavirs 🔬 ![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/OpenOmics/metavirs?color=blue&include_prereleases) [![tests](https://github.com/OpenOmics/metavirs/workflows/tests/badge.svg)](https://github.com/OpenOmics/metavirs/actions/workflows/main.yaml) [![docs](https://github.com/OpenOmics/metavirs/workflows/docs/badge.svg)](https://github.com/OpenOmics/metavirs/actions/workflows/docs.yml) [![GitHub issues](https://img.shields.io/github/issues/OpenOmics/metavirs?color=brightgreen)](https://github.com/OpenOmics/metavirs/issues)  [![GitHub license](https://img.shields.io/github/license/OpenOmics/metavirs)](https://github.com/OpenOmics/metavirs/blob/main/LICENSE) 

> **_Metagenomics Viral Sequencing Pipeline_**. This is the home of the pipeline, metavirs. Its long-term goals: to assemble, annotate, and screen enviromental samples like no pipeline before!

## Overview
Welcome to metavirs! Before getting started, we highly recommend reading through [metavirs' documentation](https://openomics.github.io/metavirs/).

The **`./metavirs`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

 * [<code>metavirs <b>run</b></code>](https://openomics.github.io/metavirs/usage/run/): Run the metavirs pipeline with your input files.
 * [<code>metavirs <b>unlock</b></code>](https://openomics.github.io/metavirs/usage/unlock/): Unlocks a previous runs output directory.
 * [<code>metavirs <b>cache</b></code>](https://openomics.github.io/metavirs/usage/cache/): Cache remote resources locally, coming soon!

metavirs is a comprehensive viral metagenomics pipeline assemble, annotate, and classify enviromental samples. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of FastQ files and can be run locally on a compute instance, on-premise using a cluster, or on the cloud (feature coming soon!). A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM, or run on AWS using Tibanna (feature coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](https://openomics.github.io/metavirs/usage/run/) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](https://openomics.github.io/metavirs/faq/questions/) prior to [opening an issue on Github](https://github.com/OpenOmics/metavirs/issues).

## Dependencies
**Requires:** `singularity>=3.5`  `snakemake>=6.0` `mamba`

At the current moment, the pipeline uses a mixture of enviroment modules, docker images, and conda environments; however, this will be changing soon! In the very near future, the pipeline will only use docker images. With that being said, [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), [singularity](https://singularity.lbl.gov/all-releases), and [mamba](https://github.com/mamba-org/mamba#installation) must be installed on the target system. Snakemake orchestrates the execution of each step in the pipeline. To guarantee the highest level of reproducibility, each step of the pipeline will rely on versioned images from [DockerHub](https://hub.docker.com/orgs/nciccbr/repositories). Snakemake uses singularity to pull these images onto the local filesystem prior to job execution, and as so, snakemake and singularity will be the only two dependencies in the future.

## Installation
Please clone this repository to your local filesystem using the following command:
```bash
# Clone Repository from Github
git clone https://github.com/OpenOmics/metavirs.git
# Change your working directory
cd metavirs/
# Add dependencies to $PATH
# Biowulf users should run
source /data/CCBR_Pipeliner/db/PipeDB/Conda/etc/profile.d/mamba.sh 
module load snakemake 
module load singularity
# Get usage information
./metvirs -h
```

## Contribute 
This site is a living document, created for and by members like you. metavirs is maintained by the members of NCBR and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository](https://github.com/OpenOmics/metavirs).


## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
