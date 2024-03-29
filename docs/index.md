<div align="center">

  <h1 style="font-size: 250%">metavirs 🔬</h1>

  <b><i>Viral Metagenomics Pipeline</i></b><br>
  <a href="https://doi.org/10.5281/zenodo.7120936">
      <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7120936.svg" alt="DOI">
  </a>
  <a href="https://github.com/OpenOmics/metavirs/releases">
    <img alt="GitHub release" src="https://img.shields.io/github/v/release/OpenOmics/metavirs?color=blue&include_prereleases">
  </a>
  <a href="https://hub.docker.com/repository/docker/skchronicles/metavirs">
    <img alt="Docker Pulls" src="https://img.shields.io/docker/pulls/skchronicles/metavirs">
  </a><br>
  <a href="https://github.com/OpenOmics/metavirs/actions/workflows/main.yaml">
    <img alt="tests" src="https://github.com/OpenOmics/metavirs/workflows/tests/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/metavirs/actions/workflows/docs.yml">
    <img alt="docs" src="https://github.com/OpenOmics/metavirs/workflows/docs/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/metavirs/issues">
    <img alt="GitHub issues" src="https://img.shields.io/github/issues/OpenOmics/metavirs?color=brightgreen">
  </a>
  <a href="https://github.com/OpenOmics/metavirs/blob/main/LICENSE">
    <img alt="GitHub license" src="https://img.shields.io/github/license/OpenOmics/metavirs">
  </a>

  <p>
    This is the home of the pipeline, metavirs. Its long-term goals: to assemble, annotate, and classify enviromental samples like no pipeline before!
  </p>

</div>  


## Overview
Welcome to our documentation! This guide is the main source of documentation for users that are getting started with the [viral metagenomics pipeline](https://github.com/OpenOmics/metavirs/). 

The **`./metavirs`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

 * [<code>metavirs <b>run</b></code>](usage/run.md): Run the metavirs pipeline with your input files.
 * [<code>metavirs <b>unlock</b></code>](usage/unlock.md): Unlocks a previous runs output directory.
 * [<code>metavirs <b>install</b></code>](usage/install.md): Download reference files locally.
 * [<code>metavirs <b>cache</b></code>](usage/cache.md): Cache software containers locally.

metavirs is a comprehensive viral metagenomics pipeline to assemble, annotate, and classify microorganisms in enviromental samples. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of FastQ files and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](usage/run.md) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](faq/questions.md) prior to [opening an issue on Github](https://github.com/OpenOmics/metavirs/issues).

## Citation

If you use this software, please cite it as below:  

=== "BibTex"

    ```
    @software{Kuhn_OpenOmics_metavirs_2022,
        author = {Kuhn, Skyler and Schaughency, Paul},
        doi = {10.5281/zenodo.7120936},
        month = {9},
        title = {{OpenOmics/metavirs}},
        url = {https://github.com/OpenOmics/metavirs/},
        year = {2022}
    }
    ```

=== "APA"

    ``` 
    Kuhn, S., & Schaughency, P. (2022). OpenOmics/metavirs [Computer software]. https://doi.org/10.5281/zenodo.7120936
    ```

For more citation style options, please visit the pipeline's [Zenodo page](https://doi.org/10.5281/zenodo.7120936).

## Contribute 

This site is a living document, created for and by members like you. metavirs is maintained by the members of OpenOmics and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository :octicons-heart-fill-24:{ .heart }](https://github.com/OpenOmics/metavirs).

## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
