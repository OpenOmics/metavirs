# <code>metavirs <b>run</b></code>

## 1. About 
The `metavirs` executable is composed of several inter-related sub commands. Please see `metavirs -h` for all available options.

This part of the documentation describes options and concepts for <code>metavirs <b>run</b></code> sub command in more detail. With minimal configuration, the **`run`** sub command enables you to start running metavirs pipeline. 

Setting up the metavirs pipeline is fast and easy! In its most basic form, <code>metavirs <b>run</b></code> only has *two required inputs*.

## 2. Synopsis
```text
$ metavirs run [--help] [--mode <slurm,local>] \
     [--job-name JOB_NAME] [--dry-run] [--silent] \
     [--singularity-cache SINGULARITY_CACHE] \
     [--sif-cache SIF_CACHE] \
     [--tmpdir TMP_DIR] \
     [--threads THREADS] \
      --input INPUT [INPUT ...] \
      --output OUTPUT
```

The synopsis for each command shows its arguments and their usage. Optional arguments are shown in square brackets.

A user **must** provide a list of FastQ (globbing is supported) to analyze via `--input` argument and an output directory to store results via `--output` argument.

Use you can always use the `-h` option for information on a specific command. 

### 2.1 Required Arguments

Each of the following arguments are required. Failure to provide a required argument will result in a non-zero exit-code.

  `--input INPUT [INPUT ...]`  
> **Input FastQ or BAM file(s).**  
> *type: file(s)*  
> 
> One or more FastQ files can be provided. The pipeline does NOT support single-end data. From the command-line, each input file should seperated by a space. Globbing is supported! This makes selecting FastQ files easy. Input FastQ files should always be gzipp-ed.
> 
> ***Example:*** `--input .tests/*.R?.fastq.gz`

---  
  `--output OUTPUT`
> **Path to an output directory.**   
> *type: path*
>   
> This location is where the pipeline will create all of its output files, also known as the pipeline's working directory. If the provided output directory does not exist, it will be created automatically.
> 
> ***Example:*** `--output /data/$USER/metavirs_out`

### 2.2 Options

Each of the following arguments are optional, and do not need to be provided. 

  `-h, --help`            
> **Display Help.**  
> *type: boolean flag*
> 
> Shows command's synopsis, help message, and an example command
> 
> ***Example:*** `--help`

---  
  `--dry-run`            
> **Dry run the pipeline.**  
> *type: boolean flag*
> 
> Displays what steps in the pipeline remain or will be run. Does not execute anything!
>
> ***Example:*** `--dry-run`

---  
  `--silent`            
> **Silence standard output.**  
> *type: boolean flag*
> 
> Reduces the amount of information directed to standard output when submitting master job to the job scheduler. Only the job id of the master job is returned.
>
> ***Example:*** `--silent`

---  
  `--mode {slurm,local}`  
> **Execution Method.**  
> *type: string*  
> *default: slurm*
> 
> Execution Method. Defines the mode or method of execution. Vaild mode options include: slurm or local. 
> 
> ***slurm***    
> The slurm execution method will submit jobs to the [SLURM workload manager](https://slurm.schedmd.com/). It is recommended running metavirs in this mode as execution will be significantly faster in a distributed environment. This is the default mode of execution.
>
> ***local***  
> Local executions will run serially on compute instance. This is useful for testing, debugging, or when a users does not have access to a high performance computing environment. If this option is not provided, it will default to a local execution mode. 
> 
> ***Example:*** `--mode slurm`

---  
  `--job-name JOB_NAME`  
> **Set the name of the pipeline's master job.**  
> *type: string*
> *default: pl:metavirs*
> 
> When submitting the pipeline to a job scheduler, like SLURM, this option always you to set the name of the pipeline's master job. By default, the name of the pipeline's master job is set to "pl:metavirs".
> 
> ***Example:*** `--job-name pl_id-42`

---  
  `--singularity-cache SINGULARITY_CACHE`  
> **Overrides the $SINGULARITY_CACHEDIR environment variable.**  
> *type: path*  
> *default: `--output OUTPUT/.singularity`*
>
> Singularity will cache image layers pulled from remote registries. This ultimately speeds up the process of pull an image from DockerHub if an image layer already exists in the singularity cache directory. By default, the cache is set to the value provided to the `--output` argument. Please note that this cache cannot be shared across users. Singularity strictly enforces you own the cache directory and will return a non-zero exit code if you do not own the cache directory! See the `--sif-cache` option to create a shareable resource. 
> 
> ***Example:*** `--singularity-cache /data/$USER/.singularity`

---  
  `--sif-cache SIF_CACHE`
> **Path where a local cache of SIFs are stored.**  
> *type: path*  
>
> Uses a local cache of SIFs on the filesystem. This SIF cache can be shared across users if permissions are set correctly. If a SIF does not exist in the SIF cache, the image will be pulled from Dockerhub and a warning message will be displayed. The `metavirs cache` subcommand can be used to create a local SIF cache. Please see `metavirs cache` for more information. This command is extremely useful for avoiding DockerHub pull rate limits. It also remove any potential errors that could occur due to network issues or DockerHub being temporarily unavailable. We recommend running metavirs with this option when ever possible.
> 
> ***Example:*** `--singularity-cache /data/$USER/SIFs`

---  
  `--threads THREADS`   
> **Max number of threads for each process.**  
> *type: int*  
> *default: 2*
> 
> Max number of threads for each process. This option is more applicable when running the pipeline with `--mode local`.  It is recommended setting this vaule to the maximum number of CPUs available on the host machine.
> 
> ***Example:*** `--threads 12`

## 3. Example
```bash 
# Step 1.) Grab an interactive node
# Do not run on head node!
sinteractive --mem=8g --cpus-per-task=4
module purge
/data/CCBR_Pipeliner/db/PipeDB/Conda/bin/conda activate base
module load singularity snakemake

# Step 2A.) Dry-run the pipeline
./metavirs run --input .tests/*.gz \
    --output /data/$USER/metavirs_out \
    --mode slurm \
    --dry-run

# Step 2B.) Run the viral metagenomics pipeline
# The slurm mode will submit jobs to the cluster.
# It is recommended running metavirs in this mode.
./metavirs run --input .tests/*.gz \
    --output /data/$USER/metavirs_out \
    --mode slurm 
```