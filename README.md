<div align="center">
   
  <h1>metavirs 🔬</h1>
  
  **_Viral Metagenomics Pipeline_**

  [![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/OpenOmics/metavirs?color=blue&include_prereleases)](https://github.com/OpenOmics/metavirs/releases) [![Docker Pulls](https://img.shields.io/docker/pulls/skchronicles/metavirs)](https://hub.docker.com/repository/docker/skchronicles/metavirs) [![tests](https://github.com/OpenOmics/metavirs/workflows/tests/badge.svg)](https://github.com/OpenOmics/metavirs/actions/workflows/main.yaml) [![docs](https://github.com/OpenOmics/metavirs/workflows/docs/badge.svg)](https://github.com/OpenOmics/metavirs/actions/workflows/docs.yml) [![GitHub issues](https://img.shields.io/github/issues/OpenOmics/metavirs?color=brightgreen)](https://github.com/OpenOmics/metavirs/issues)  [![GitHub license](https://img.shields.io/github/license/OpenOmics/metavirs)](https://github.com/OpenOmics/metavirs/blob/main/LICENSE)  
  
  <i>
    This is the home of the pipeline, metavirs. Its long-term goals: to assemble, annotate, and screen enviromental samples like no pipeline before!
  </i>
</div>

## Overview
Welcome to metavirs! Before getting started, we highly recommend reading through [metavirs' documentation](https://openomics.github.io/metavirs/).

The **`./metavirs`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

 * [<code>metavirs <b>run</b></code>](https://openomics.github.io/metavirs/usage/run/): Run the metavirs pipeline with your input files.
 * [<code>metavirs <b>unlock</b></code>](https://openomics.github.io/metavirs/usage/unlock/): Unlocks a previous runs output directory.
 * [<code>metavirs <b>install</b></code>](https://openomics.github.io/metavirs/usage/install/): Download reference files locally.
 * [<code>metavirs <b>cache</b></code>](https://openomics.github.io/metavirs/usage/cache/): Cache software containers locally.


metavirs is a comprehensive viral metagenomics pipeline to assemble, annotate, and classify microorganisms in enviromental samples. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of FastQ files and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](https://openomics.github.io/metavirs/usage/run/) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](https://openomics.github.io/metavirs/faq/questions/) prior to [opening an issue on Github](https://github.com/OpenOmics/metavirs/issues).

## Dependencies
**Requires:** `singularity>=3.5`  `snakemake>=6.0` 

At the current moment, the pipeline only has two dependencies: snakemake and singularity. With that being said, [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and [singularity](https://singularity.lbl.gov/all-releases) must be installed on the target system. Snakemake orchestrates the execution of each step in the pipeline. To guarantee the highest level of reproducibility, each step of the pipeline relies on versioned images from [DockerHub](https://hub.docker.com/repository/docker/skchronicles/metavirs). 

## Installation
Please clone this repository to your local filesystem using the following command:
```bash
# Clone Repository from Github
git clone https://github.com/OpenOmics/metavirs.git
# Change your working directory
cd metavirs/
# Add dependencies to $PATH
# Biowulf users should run
module load snakemake singularity
# Get usage information
./metvirs -h
```

## Contribute 
This site is a living document, created for and by members like you. metavirs is maintained by the members of OpenOmics and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository](https://github.com/OpenOmics/metavirs).


## Citation

If you use this software, please cite it as below:  

**BibTex**  
```text
@software{Kuhn_OpenOmics_metavirs_2022,
    author = {Kuhn, Skyler and Schaughency, Paul},
    doi = {10.5281/zenodo.7120936},
    month = {9},
    title = {{OpenOmics/metavirs}},
    url = {https://github.com/OpenOmics/metavirs/},
    year = {2022}
}
```

**APA**  
```text
Kuhn, S., & Schaughency, P. (2022). OpenOmics/metavirs [Computer software]. https://doi.org/10.5281/zenodo.7120936
```

For more citation style options, please visit the pipeline's [Zenodo page](https://doi.org/10.5281/zenodo.7120936).

## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
