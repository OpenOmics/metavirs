# Paired-end snakemake rules imported in the main Snakefile.
from scripts.common import abstract_location, references, allocated


# Pre alignment QC-related rules
rule validator:
    """
    Quality-control step to ensure the input FastQC files are not corrupted or
    incomplete prior to running the entire workflow. Please note this rule will
    only run if the --use-singularity flag is provided to snakemake.
    @Input:
        Raw FastQ file (scatter)
    @Output:
        Log file containing any warnings or errors on file
    """
    input:
        r1=join(workpath,"{name}.R1.fastq.gz"),
        r2=join(workpath,"{name}.R2.fastq.gz"),
    output:
        log1=join(workpath,"rawQC","{name}.validated.R1.fastq.log"),
        log2=join(workpath,"rawQC","{name}.validated.R2.fastq.log"),
    params:
        rname='validfq',
        outdir=join(workpath,"rawQC"),
    container: config['images']['fastqvalidator']
    shell: """
    mkdir -p {params.outdir}
    fastQValidator --noeof \\
        --file {input.r1} > {output.log1}
    fastQValidator --noeof \\
        --file {input.r2} > {output.log2}
    """


rule rawfastqc:
    """
    Quality-control step to assess sequencing quality of the raw data prior removing
    adapter sequences. FastQC generates a set of basic statistics to identify problems
    that can arise during sequencing or library preparation.
    @Input:
        List of Raw FastQ files (gather)
    @Output:
        List of FastQC reports and zip file containing data quality information
    """
    input:
        expand(join(workpath,"{name}.R1.fastq.gz"), name=samples),
        expand(join(workpath,"{name}.R2.fastq.gz"), name=samples),
    output:
        expand(join(workpath,"rawQC","{name}.R1_fastqc.zip"), name=samples),
        expand(join(workpath,"rawQC","{name}.R2_fastqc.zip"), name=samples),
    params:
        rname='rawfqc',
        outdir=join(workpath,"rawQC"),
    threads: int(allocated("threads", "rawfastqc", cluster))
    envmodules: config['tools']['fastqc']
    container: config['images']['fastqc']
    shell: """
    fastqc {input} \\
        -t {threads} \\
        -o {params.outdir}
    """


# Data processing rules
rule trim:
    """
    Data-processing step to remove adapter sequences and perform quality trimming
    prior to alignment the reference genome.  Adapters are composed of synthetic
    sequences and should be removed prior to alignment.
    @Input:
        Raw FastQ file (scatter)
    @Output:
        Trimmed FastQ file
    """
    input:
        r1=join(workpath,"{name}.R1.fastq.gz"),
        r2=join(workpath,"{name}.R2.fastq.gz"),
    output:
        r1=join(workpath,"trim","{name}.R1.trim.fastq"),
        r2=join(workpath,"trim","{name}.R2.trim.fastq"),
    params:
        rname='trimfq',
        adapters=config['references']['adapters'],
    threads: int(allocated("threads", "trim", cluster))
    envmodules: config['tools']['cutadapt']
    container: config['images']['cutadapt']
    shell: """
    cutadapt -j {threads} \\
        --pair-filter=any \\
        --nextseq-trim=2 \\
        --trim-n \\
        -n 5 -O 5 \\
        -q 10,10 \\
        -m 35:35 \\
        -b file:{params.adapters} \\
        -B file:{params.adapters} \\
        -o {output.r1} -p {output.r2} \\
        {input.r1} {input.r2}
    """


rule fastqc:
    """
    Quality-control step to assess sequencing quality of the raw data after removing
    adapter sequences. This step is run after trim_pe rule. FastQC is run after adapter
    trimming to evalute if the adapter sequences were properly removed.
    @Input:
        List of Trimmed FastQ files (gather)
    @Output:
        List of FastQC reports and zip file containing data quality information
    """
    input:
        expand(join(workpath,"trim","{name}.R1.trim.fastq"), name=samples),
        expand(join(workpath,"trim","{name}.R2.trim.fastq"), name=samples),
    output:
        expand(join(workpath,"QC","{name}.R1.trim_fastqc.zip"), name=samples),
        expand(join(workpath,"QC","{name}.R2.trim_fastqc.zip"), name=samples),
    params:
        rname='fqc',
        outdir=join(workpath,"QC"),
    threads: int(allocated("threads", "fastqc", cluster))
    envmodules: config['tools']['fastqc']
    container: config['images']['fastqc']
    shell: """
    fastqc {input} \\
        -t {threads} \\
        -o {params.outdir}
    """
