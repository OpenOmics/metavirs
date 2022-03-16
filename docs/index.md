# metavirs 🔬  [![docs](https://github.com/OpenOmics/metavirs/workflows/docs/badge.svg)](https://github.com/OpenOmics/metavirs/actions) [![GitHub issues](https://img.shields.io/github/issues/OpenOmics/metavirs?color=brightgreen)](https://github.com/OpenOmics/metavirs/issues)  [![GitHub license](https://img.shields.io/github/license/OpenOmics/metavirs)](https://github.com/OpenOmics/metavirs/blob/main/LICENSE) 

> **_Metagenomics Viral Sequencing Pipeline_**. This is the home of the pipeline, metavirs. Its long-term goals: to assemble, annotate, and classify enviromental samples like no pipeline before!

---
## Overview
Welcome to metavirs' documentation! This guide is the main source of documentation for users that are getting started with the [Viral metagenomics pipeline](https://github.com/OpenOmics/metavirs/). 

The **`./metavirs`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

 * [<code>metavirs <b>run</b></code>](usage/run.md): Run the GATK4 WES pipeline with your input files.
 * [<code>metavirs <b>unlock</b></code>](usage/unlock.md): Unlocks a previous runs output directory.
 * [<code>metavirs <b>cache</b></code>](usage/cache.md): Cache remote resources locally, coming soon!

metavirs is a comprehensive viral metagenomics pipeline assemble, annotate, and classify enviromental samples. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of FastQ files and can be run locally on a compute instance, on-premise using a cluster, or on the cloud (feature coming soon!). A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM, or run on AWS using Tibanna (feature coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](usage/run.md) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](faq/questions.md) prior to [opening an issue on Github](https://github.com/OpenOmics/metavirs/issues).

## Contribute 

This site is a living document, created for and by members like you. metavirs is maintained by the members of NCBR and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository :octicons-heart-fill-24:{ .heart }](https://github.com/OpenOmics/metavirs).


## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
