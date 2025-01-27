# Paired-end snakemake rules imported in the main Snakefile.
from scripts.common import abstract_location, references, allocated

# Helper functions for running rules with
# mixed input for single-end and paired-end
# datasets 
def get_r2_raw_fastq(wildcards):
    """
    Returns a samples R2 fastq file
    if it is paired-end data.
    """
    r2 = nends[wildcards.name]
    if r2:
        # Runs in paired-end mode
        return join(workpath,"{0}.R2.fastq.gz".format(r2))
    else:
        # Runs in single-end mode
        return []


def get_r2_trim_fastq(wildcards):
    """
    Returns a samples trimmed R2 fastq file
    if it is paired-end data.
    """
    r2 = nends[wildcards.name]
    if r2:
        # Runs in paired-end mode
        return join(workpath,"{0}".format(r2),"trim","{0}.R2.trim.fastq".format(r2))
    else:
        # Runs in single-end mode
        return []


def get_r2_host_removed_fastq(wildcards):
    """
    Returns a samples trimmed, host removed R2 fastq file
    if it is paired-end data.
    """
    r2 = nends[wildcards.name]
    if r2:
        # Runs in paired-end mode
        return join(workpath,"{0}".format(r2),"trim","{0}.R2.trim.host_removed.fastq.gz".format(r2))
    else:
        # Runs in single-end mode
        return []


def get_contigs_fasta(wildcards):
    """
    Returns a samples annotated contigs fasta file.
    Metaspades does NOT support single-end data, so
    metaspades will conditionally run for a given set
    of samples. This function resolves the correct set 
    of files for a given sample based on its paired-end
    status.
    """
    # Resolves to basename of sample
    r2 = nends[wildcards.name]
    if r2:
        # Runs in paired-end mode,
        # returns metaspades and megahit
        i = [
            join(workpath,"{0}".format(wildcards.name),"output","{0}.metaspades.contigs.fa".format(wildcards.name)),
            join(workpath,"{0}".format(wildcards.name),"output","{0}.megahit.contigs.fa".format(wildcards.name))
        ]
        return i
    else:
        # Runs in single-end mode,
        # does not add metaspades to inputs list
        i = [
            join(workpath,"{0}".format(wildcards.name),"output","{0}.megahit.contigs.fa".format(wildcards.name))
        ]
        return i


def get_annotated_contigs(wildcards):
    """
    Returns a samples annotated contigs text files.
    Metaspades does NOT support single-end data, so
    metaspades will conditionally run for a given set
    of samples. This function resolves the correct set 
    of files for a given sample based on its paired-end
    status.
    """
    # Resolves to basename of sample
    r2 = nends[wildcards.name]
    if r2:
        # Runs in paired-end mode,
        # returns all assembler/annotator outputs
        i = [
            join(workpath,"{0}".format(wildcards.name),"temp","{0}.metaspades_CAT.txt".format(wildcards.name)),
            join(workpath,"{0}".format(wildcards.name),"temp","{0}.metaspades_kraken2.txt".format(wildcards.name)),
            join(workpath,"{0}".format(wildcards.name),"temp","{0}.metaspades_kraken2_filtered.txt".format(wildcards.name)),
            join(workpath,"{0}".format(wildcards.name),"temp","{0}.megahit_CAT.txt".format(wildcards.name)),
            join(workpath,"{0}".format(wildcards.name),"temp","{0}.megahit_kraken2.txt".format(wildcards.name)),
            join(workpath,"{0}".format(wildcards.name),"temp","{0}.megahit_kraken2_filtered.txt".format(wildcards.name))
        ]
        return i
    else:
        # Runs in single-end mode,
        # does not add metaspades to inputs list
        wildcards.name
        i = [
            join(workpath,"{0}".format(wildcards.name),"temp","{0}.megahit_CAT.txt".format(wildcards.name)),
            join(workpath,"{0}".format(wildcards.name),"temp","{0}.megahit_kraken2.txt".format(wildcards.name)),
            join(workpath,"{0}".format(wildcards.name),"temp","{0}.megahit_kraken2_filtered.txt".format(wildcards.name))
        ]
        return i


def get_aggregated_results(wildcards):
    """
    Returns all samples aggregated results.
    Metaspades does NOT support single-end data, so
    metaspades will conditionally run for a given set
    of samples. This function resolves the correct set 
    of files for a given sample based on its paired-end
    status.
    """
    if pe_samples:
        # Return metaspades and megahit results
        return [
            join(workpath,"Project","temp","metaspades_kraken2.krona"),
            join(workpath,"Project","temp","metaspades_kraken2_filtered.krona"),
            join(workpath,"Project","temp","metaspades_CAT.krona"),
            join(workpath,"Project","temp","megahit_kraken2.krona"),
            join(workpath,"Project","temp","megahit_kraken2_filtered.krona"),
            join(workpath,"Project","temp","megahit_CAT.krona")
        ]
    else:
        # All samples are single-end,
        # return only megahit results 
        return [
            join(workpath,"Project","temp","megahit_kraken2.krona"),
            join(workpath,"Project","temp","megahit_kraken2_filtered.krona"),
            join(workpath,"Project","temp","megahit_CAT.krona")
        ]


# Pre alignment QC-related rules
rule validator_r1:
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
    output:
        log=join(workpath,"{name}","rawQC","{name}.validated.R1.fastq.log"),
    params:
        rname='validfqr1',
        outdir=join(workpath,"{name}","rawQC"),
    threads: int(allocated("threads", "validator_r1", cluster))
    container: config['images']['fastqvalidator']
    shell: """
    mkdir -p {params.outdir}
    fastQValidator --noeof \\
        --file {input.r1} \\
    > {output.log}
    """


rule validator_r2:
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
        r2=join(workpath,"{name}.R2.fastq.gz"),
    output:
        log=join(workpath,"{name}","rawQC","{name}.validated.R2.fastq.log"),
    params:
        rname='validfqr2',
        outdir=join(workpath,"{name}","rawQC"),
    threads: int(allocated("threads", "validator_r2", cluster))
    container: config['images']['fastqvalidator']
    shell: """
    mkdir -p {params.outdir}
    fastQValidator --noeof \\
        --file {input.r2} \\
    > {output.log}
    """


rule rawfastqc:
    """
    Quality-control step to assess sequencing quality of the raw data prior removing
    adapter sequences. FastQC generates a set of basic statistics to identify problems
    that can arise during sequencing or library preparation.
    @Input:
        Raw FastQ file (scatter)
    @Output:
        FastQC report and zip file containing data quality information
    """
    input:
        join(workpath,"{name}.R1.fastq.gz"),
        get_r2_raw_fastq,
    output:
        join(workpath,"{name}","rawQC","{name}.R1_fastqc.zip"),
    params:
        rname='rawfqc',
        outdir=join(workpath,"{name}","rawQC"),
        tmpdir=tmpdir,
    threads: int(allocated("threads", "rawfastqc", cluster))
    # envmodules: config['tools']['fastqc']
    container: config['images']['metavirs']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT
    
    # Running fastqc with local
    # disk or a tmpdir, fastqc
    # has been observed to lock
    # up gpfs filesystems, adding
    # this on request by HPC staff. 
    fastqc {input} \\
        -t {threads} \\
        -o "${{tmp}}"
    
    # Copy output files from tmpdir
    # to output directory
    find "${{tmp}}" \\
        -type f \\
        \\( -name '*.html' -o -name '*.zip' \\) \\
        -exec cp {{}} {params.outdir} \\;
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
        r2=get_r2_raw_fastq,
    output:
        r1=join(workpath,"{name}","trim","{name}.R1.trim.fastq"),
    params:
        rname='trimfq',
        adapters=config['references']['adapters'],
        # Building single-end and paired-end options
        # --pair-filter: PE='--pair-filter=any', SE=''
        pair_filter_option = lambda w: "--pair-filter=any" \
            if nends[w.name] else "",
        # -m: PE='-m 35:35', SE='-m 35'
        m_option = lambda w: "-m {0}".format(
            "35:35"
        ) if nends[w.name] else "-m 35",
        # -B: PE='-B file:{params.adapters}', SE=''
        B_option = lambda w: "-B file:{0}".format(
            config['references']['adapters']
        ) if nends[w.name] else "",
        # -p: PE='-p {output.r2}', SE=''
        p_option = lambda w: "-p {0}.R2.trim.fastq".format(
            join(workpath, w.name, "trim", w.name)
        ) if nends[w.name] else "",
        # Input R2 FastQ file: PE='{input.r2}', SE=''
        r2_option = lambda w: "{0}.R2.fastq.gz".format(
            join(workpath, w.name)
        ) if nends[w.name] else "",
    threads: int(allocated("threads", "trim", cluster))
    # envmodules: config['tools']['cutadapt']
    container: config['images']['metavirs']
    shell: """
    cutadapt -j {threads} {params.pair_filter_option} \\
        --nextseq-trim=2 \\
        --trim-n \\
        -n 5 -O 5 \\
        -q 10,10 \\
        {params.m_option} \\
        -b file:{params.adapters} {params.B_option} \\
        -o {output.r1} {params.p_option} \\
        {input.r1} {params.r2_option}
    """


rule fastqc:
    """
    Quality-control step to assess sequencing quality of the raw data after removing
    adapter sequences. This step is run after trim_pe rule. FastQC is run after adapter
    trimming to evalute if the adapter sequences were properly removed.
    @Input:
        Trimmed FastQ file (scatter)
    @Output:
        FastQC report and zip file containing data quality information
    """
    input:
        r1=join(workpath,"{name}","trim","{name}.R1.trim.fastq"),
    output:
        join(workpath,"{name}","QC","{name}.R1.trim_fastqc.zip"),
    params:
        rname='fqc',
        outdir=join(workpath,"{name}","QC"),
        tmpdir=tmpdir,
        # Building single-end and paired-end options
        # Input trimmed R2 FastQ file: PE='{trim.r2}', SE=''
        r2_option = lambda w: "{0}.R2.trim.fastq".format(
            join(workpath, w.name, "trim", w.name)
        ) if nends[w.name] else "",
    threads: int(allocated("threads", "fastqc", cluster))
    # envmodules: config['tools']['fastqc']
    container: config['images']['metavirs']
    shell: """
    # Setups temporary directory for
    # intermediate files with built-in 
    # mechanism for deletion on exit
    if [ ! -d "{params.tmpdir}" ]; then mkdir -p "{params.tmpdir}"; fi
    tmp=$(mktemp -d -p "{params.tmpdir}")
    trap 'rm -rf "${{tmp}}"' EXIT
    
    # Running fastqc with local
    # disk or a tmpdir, fastqc
    # has been observed to lock
    # up gpfs filesystems, adding
    # this on request by HPC staff. 
    fastqc {input.r1} {params.r2_option} \\
        -t {threads} \\
        -o "${{tmp}}"
    
    # Copy output files from tmpdir
    # to output directory
    find "${{tmp}}" \\
        -type f \\
        \\( -name '*.html' -o -name '*.zip' \\) \\
        -exec cp {{}} {params.outdir} \\;
    """


rule remove_host:
    """
    Data-processing step to remove any sources of contamination from host. Reads
    are aligned against a host bowtie2 index. Any reads that align against the 
    host reference genome are then removed from the trimmed FastQ files. 
    @Input:
        Trimmed FastQ file (scatter)
    @Output:
        Trimmed, host remove FastQ file
    """
    input:
        r1=join(workpath,"{name}","trim","{name}.R1.trim.fastq"),
    output:
        sam=temp(join(workpath,"{name}","temp","{name}.host_contaminated.sam")),
        bam=temp(join(workpath,"{name}","temp","{name}.host_contaminated.bam")),
        sorted_bam=temp(join(workpath,"{name}","temp","{name}.host_contaminated.sorted.bam")),
        r1=join(workpath,"{name}","trim","{name}.R1.trim.host_removed.fastq.gz"),
    params:
        rname='rmhost',
        host_index=join(config['references']['host_bowtie2_index'], 'Hosts'),
        # Building single-end and paired-end options
        # Input trimmed R2 FastQ file: PE='-2 {trim.r2}', SE=''
        r2_option = lambda w: "-2 {0}.R2.trim.fastq".format(
            join(workpath, w.name, "trim", w.name)
        ) if nends[w.name] else "",
        # Samtools view flag: PE='-f 12', SE='-f 4'
        # https://broadinstitute.github.io/picard/explain-flags.html
        f_option = lambda w: "-f 12" if nends[w.name] else "-f 4",
        # Output host removed R2 fastq: PE='-2 output.r2', SE=''
        r2_output = lambda w: "-2 {0}.R2.trim.host_removed.fastq.gz".format(
            join(workpath, w.name, "trim", w.name)
        ) if nends[w.name] else "",
        # Bowtie2 SE vs PE conditional option: PE='-1', SE='-U'
        bowtie2_sepe = lambda w: "-1" if nends[w.name] else "-U",
        # Samtools fastq output: PE='-1 {output.r1}', SE='-0 {output.r1}'
        samtools_sepe = lambda w: "-1" if nends[w.name] else "-0",
        # Samtools fastq, PE options: PE='-0 /dev/null -s /dev/null', SE=''
        samtools_pe = lambda w: "-0 {0} -s {0}".format(
            join(os.path.sep, 'dev', 'null')
        ) if nends[w.name] else "",
    threads: int(allocated("threads", "remove_host", cluster))
    # envmodules: config['tools']['bowtie'], config['tools']['samtools']
    container: config['images']['metavirs']
    shell: """
    # Align against host to remove contamination
    bowtie2 -p {threads} \\
        -x {params.host_index} \\
        {params.bowtie2_sepe} {input.r1} {params.r2_option} \\
        -S {output.sam}
    samtools view -@ {threads} \\
        -bS {params.f_option} -F 256 \\
        {output.sam} > {output.bam}
    samtools sort -n -@ {threads} \\
        {output.bam} \\
        -o {output.sorted_bam}
    samtools fastq -@ {threads} \\
        {output.sorted_bam} \\
        {params.samtools_sepe} {output.r1} {params.r2_output} \\
        -n {params.samtools_pe}
    """


rule kraken_viral:
    """
    Data-processing step for taxonomic classification of microbes. Kraken2
    uses exact k-mer matches to achieve high accuracy and fast classification
    speeds. This classifier matches each k-mer within a query sequence to the
    lowest common ancestor (LCA) of all genomes containing the given k-mer. 
    @Input:
        Trimmed, host remove FastQ file (scatter)
    @Output:
        Taxonomic classification of trimmed, host remove reads
    """
    input:
        r1=join(workpath,"{name}","trim","{name}.R1.trim.host_removed.fastq.gz"),
    output:
        report=join(workpath,"{name}","info","{name}.reads_kraken2_report.txt"),
        k2txt=join(workpath,"{name}","kraken2","{name}.reads.kraken2"),
        kronatxt=join(workpath,"{name}","kraken2","{name}.reads.kraken2.krona"),
        html=join(workpath,"{name}","kraken2","{name}.reads.krona.html"),
    params:
        rname='krakenviral',
        viral_db=config['references']['kraken2_viral_db'],
        krona_ref=config['references']['kronatools'],
        # Building single-end and paired-end options
        # Kraken PE flag: PE='--paired', SE=''
        paired_option = lambda w: "--paired" if nends[w.name] else "",
        # Input host removed R2 fastq: PE='input.r2', SE=''
        r2_option = lambda w: "{0}.R2.trim.host_removed.fastq.gz".format(
            join(workpath, w.name, "trim", w.name)
        ) if nends[w.name] else "",
    threads: int(allocated("threads", "kraken_viral", cluster))
    # envmodules: config['tools']['kraken'], config['tools']['kronatools']
    container: config['images']['metavirs']
    shell: """
    # Run kraken against viral database
    kraken2 --threads {threads} \\
        --db {params.viral_db} \\
        {params.paired_option} {input.r1} {params.r2_option} \\
        --report {output.report} \\
        > {output.k2txt}
    awk -v OFS='\\t' '{{if ($1 == "C") print $2,$3}}' {output.k2txt} \\
        > {output.kronatxt}
    ktImportTaxonomy \\
        -tax {params.krona_ref} \\
        -o {output.html} \\
        {output.kronatxt}
    """


rule metaspades:
    """
    Data-processing step to assembly reads for metagenomic analysis. 
    @Input:
        Trimmed, host remove FastQ file (scatter)
    @Output:
        Kraken report and annotation report of assembled contigs,
        and aligned reads against the assembled contigs. 
    """
    input:
        r1=join(workpath,"{name}","trim","{name}.R1.trim.host_removed.fastq.gz"),
    output:
        contigs=join(workpath,"{name}","output","{name}.metaspades.contigs.fa"),
        report=join(workpath,"{name}","info","{name}.metaspades.contigs_kraken2_report.txt"),
        k2txt=join(workpath,"{name}","kraken2","{name}.metaspades.contigs.kraken2"),
        krona=join(workpath,"{name}","kraken2","{name}.metaspades.contigs.kraken2.krona"),
        tmp1=join(workpath,"{name}","temp","{name}.metaspades_kraken2.txt"),
        ktaxids=join(workpath,"{name}","temp","{name}.metaspades_kraken2_taxid.txt"),
        ktmp=join(workpath,"{name}","temp","{name}.metaspades_kraken2.temp"),
        knames=join(workpath,"{name}","temp","{name}.metaspades_kraken2_names.txt"),
        kcounts_temp=join(workpath,"{name}","temp","{name}.metaspades_kraken2_counts.temp"),
        kcounts_txt=join(workpath,"{name}","temp","{name}.metaspades_kraken2_counts.txt"),
        kflt=join(workpath,"{name}","temp","{name}.metaspades_kraken2_filtered.txt"),
        ktaxlineage=join(workpath,"{name}","temp","{name}.metaspades_kraken2_taxlineage.tsv"),
        kviral_unf=join(workpath,"{name}","output","{name}.metaspades_kraken2_unfiltered_viraltable.tsv"),
        kviral_flt=join(workpath,"{name}","output","{name}.metaspades_kraken2_filtered_viraltable.tsv"),
        kviral_flt_family=join(workpath,"{name}","output","{name}.metaspades_kraken2_filtered_viraltable.family-level.tsv"),
        kviral_flt_genus=join(workpath,"{name}","output","{name}.metaspades_kraken2_filtered_viraltable.genus-level.tsv"),
        kviral_flt_species=join(workpath,"{name}","output","{name}.metaspades_kraken2_filtered_viraltable.species-level.tsv"),
        cat_class=join(workpath,"{name}","CAT","{name}.metaspades.contig2classification.txt"),
        cat_names=join(workpath,"{name}","CAT","{name}.metaspades.official_names.txt"),
        cat_summary=join(workpath,"{name}","CAT","{name}.metaspades.summary.txt"),
        taxids=join(workpath,"{name}","CAT","{name}.metaspades.contig_taxids.txt"),
        tmp2=join(workpath,"{name}","temp","{name}.metaspades_CAT.txt"),
        tmp3=join(workpath,"{name}","temp","{name}.metaspades.kraken2.viral.names.txt"),
        kraken_contigs=join(workpath,"{name}","output","{name}.metaspades.kraken2_viral.contigs.fa"),
        tmp4=join(workpath,"{name}","temp","{name}.metaspades.CAT.viral.names.txt"),
        cat_contigs=join(workpath,"{name}","output","{name}.metaspades.cat_viral.contigs.fa"),
        sam=temp(join(workpath,"{name}","temp","{name}.metaspades.sam")),
        bam=temp(join(workpath,"{name}","temp","{name}.metaspades.bam")),
        final=join(workpath,"{name}","output","{name}.metaspades.bam"),        
    params:
        rname='metaspades',
        script=join(workpath, "workflow", "scripts", "collapse_on_taxa.py"),
        filter_length=config['options']['length_filter'],
        taxonkit_db=config['references']['taxonkit_db'],
        viral_db=config['references']['kraken2_viral_db'],
        cat_db=config['references']['CAT_db'],
        cat_tax=config['references']['CAT_taxonomy'],
        cat_dep=config['references']['CAT_diamond'],
        cat_dir=join(workpath,"{name}","CAT"),
        memory=allocated("mem", "metaspades", cluster).lower().rstrip('g'),
        # Building single-end and paired-end options
        # Metaspade SE vs PE conditional option: PE='-1', SE='-s'
        metaspades_sepe = lambda w: "-1" if nends[w.name] else "-s",
        # Input host removed R2 fastq: PE='-2 input.r2', SE=''
        r2_option = lambda w: "-2 {0}.R2.trim.host_removed.fastq.gz".format(
            join(workpath, w.name, "trim", w.name)
        ) if nends[w.name] else "",
        # Bowtie2 SE vs PE conditional option: PE='-1', SE='-U'
        bowtie2_sepe = lambda w: "-1" if nends[w.name] else "-U",
    threads: int(allocated("threads", "metaspades", cluster))
    # conda: config['conda']['CAT']
    # envmodules: 
    #    config['tools']['spades'], config['tools']['kraken'], config['tools']['kronatools'],
    #    config['tools']['seqtk'], config['tools']['bowtie'], config['tools']['samtools']
    container: config['images']['metavirs']
    shell: """
    set -x
    mkdir -p {wildcards.name}/metaspades
    # metaspades output directory cannot
    # exist prior to running
    if [ -d "{wildcards.name}/metaspades" ]; then
        rm -rf "{wildcards.name}/metaspades"
    fi

    # Assemble reads into contigs with metaspades
    metaspades.py -t {threads} \\
        -m {params.memory} \\
        {params.metaspades_sepe} {input.r1} {params.r2_option} \\
        -o {wildcards.name}/metaspades
    mkdir -p {wildcards.name}/output
    mv {wildcards.name}/metaspades/contigs.fasta {output.contigs}

    # Annotate the contigs with kraken2
    kraken2  --threads {threads} \\
        --db {params.viral_db} {output.contigs} \\
        --report {output.report} \\
    > {output.k2txt}
    awk -v OFS='\\t' '{{split($2,a,"_"); if ($1 == "C") print $2,$3,a[6]}}' \\
        {output.k2txt} > {output.krona}
    cp {output.krona} {output.tmp1}

    # Check the number of assembled contigs
    # prior to running CAT. CAT contigs can
    # fail if this file is empty or only 
    # contains a few contigs.
    mkdir -p {params.cat_dir}
    n_contigs=$(grep -c '^>' {output.contigs} || true)
    echo "Number of contigs from metaspades: ${{n_contigs}}"
    {{ CAT contigs -n {threads} \\
        --force \\
        -c {output.contigs} \\
        -d {params.cat_db} \\
        -t {params.cat_tax} \\
        --out_prefix {params.cat_dir}/{wildcards.name}.metaspades \\
        --path_to_diamond {params.cat_dep}
    CAT add_names -i {output.cat_class} \\
        -o {output.cat_names} \\
        -t {params.cat_tax} \\
        --only_official
    CAT summarise \\
        -c {output.contigs} \\
        -i {output.cat_names} \\
        -o {output.cat_summary}
    }} || {{
        # CAT can fail if provided an empty FASTA file
        echo "WARNING: CAT failed!" 
        echo "This could be due to an issue upstream,"
        echo "such as insufficent sequencing depth,"
        echo "low viral load, or poor quality reads."
        echo "Each of these issues can lead to poor"
        echo "quality assemblies and a low number of"
        echo "assembled contigs."
        echo "Please check the input and output of"
        echo "metaspades to troubleshoot the issue!"
        touch {output}
    }}

    # Try to grep for Virsuses, may pipefail
    # if nothing is found, touch output file 
    {{ grep "Viruses:" {output.cat_names} \\
        | awk -v OFS='\\t' '{{split($8,a,";"); \\
            split(a[split($8,a,";")],b,"*"); \\
            split($1,c,"_"); print $1,b[1],c[6]}}' \\
        > {output.taxids}
    }} || touch {output.taxids} 
    cp {output.taxids} {output.tmp2}

    cut -f1 {output.krona} > {output.tmp3}
    seqtk subseq \\
        {output.contigs} \\
        {output.tmp3} \\
    > {output.kraken_contigs}
    cut -f1 {output.taxids} > {output.tmp4}
    seqtk subseq \\
        {output.contigs} \\
        {output.tmp4} \\
    > {output.cat_contigs}

    # Build an index with assembled contigs
    # and align reads against them
    {{ bowtie2-build --threads {threads} \\
        {output.contigs} \\
        {wildcards.name}/temp/{wildcards.name}_metaspades
    bowtie2 -p {threads} \\
        -x {wildcards.name}/temp/{wildcards.name}_metaspades \\
        {params.bowtie2_sepe} {input.r1} {params.r2_option} \\
        -S {output.sam}
    samtools view -@ {threads} \\
        -bS {output.sam} > {output.bam} 
    samtools sort \\
        -@ {threads} \\
        {output.bam}  \\
        -o {output.final}
    samtools index \\
        -@ {threads} {output.final}
    }} || {{
        # Bowtie failed, maybe due to empty contigs FASTA file
        echo "WARNING: An error occurred while running bowtie!"
        echo "Please see the check the inputs to bowtie2-build"
        echo "and bowtie2 to debug the issue."
        touch {output}
    }}
    
    # Filter results based on contig length,
    # get counts for the number of reads aligning
    # to each contig, and create a counts table
    # with counts and taxonomic classifications,
    # Get taxids from intermediate results
    cut -f2 {output.tmp1} \\
    > {output.ktaxids}
    # Get contig names and coverage
    cut -f1,3 {output.tmp1} \\
    > {output.ktmp}
    # Get contig names
    cut -f1 {output.tmp1} \\
    > {output.knames}
    # Get counts of each contig in matching order
    set +o pipefail  # grep will fail -f file is empty
    samtools idxstats {output.final} \\
        | grep -f {output.knames} \\
        | cut -f3 \\
    > {output.kcounts_temp}
    set -o pipefail
    # Create text file with contig name, 
    # coverage, and counts
    paste \\
        {output.ktmp} \\
        {output.kcounts_temp} \\
    > {output.kcounts_txt}
    # Create a unfiltered and contig length 
    # filtered table with taxonomic information
    awk -F '_' -v OFS='\\t' '$4+0>={params.filter_length} {{print}}' \\
        {output.tmp1} \\
    > {output.kflt}
    cat {output.ktaxids} \\
        | taxonkit reformat \\
            --data-dir {params.taxonkit_db} \\
            -I 1 \\
            -r "Unassigned" \\
            -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}\\t{{t}}" \\
    > {output.ktaxlineage}
    # Unfiltered viral taxonomic table
    echo -e "contig\\tcov\\tcount\\ttaxid\\tkingdom\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\tstrain" \\
    > {output.kviral_unf}
    paste \\
        {output.kcounts_txt} \\
        {output.ktaxlineage} \\
    >> {output.kviral_unf}
    # Filtered viral taxonomic table
    head -1 {output.kviral_unf} \\
    > {output.kviral_flt}
    awk -F '_' -v OFS='\\t' \\
        '{{ if (NR==1) {{print}} else if ($4+0>={params.filter_length}) {{print}} }}' \\
        {output.kviral_unf} \\
    >> {output.kviral_flt}
    # Collapse and aggregate cov/counts at family-level
    python3 {params.script} \\
        --input {output.kviral_flt} \\
        --output {output.kviral_flt_family} \\
        --collapse-on-taxa family
    # Collapse and aggregate cov/counts at genus-level
    python3 {params.script} \\
        --input {output.kviral_flt} \\
        --output {output.kviral_flt_genus} \\
        --collapse-on-taxa genus
    # Collapse and aggregate cov/counts at species-level
    python3 {params.script} \\
        --input {output.kviral_flt} \\
        --output {output.kviral_flt_species} \\
        --collapse-on-taxa species
    """


rule megahit:
    """
    Data-processing step to assembly reads for metagenomic analysis. 
    @Input:
        Trimmed, host remove FastQ file (scatter)
    @Output:
        Kraken report and annotation report of assembled contigs,
        and aligned reads against the assembled megahit contigs. 
    """
    input:
        r1=join(workpath,"{name}","trim","{name}.R1.trim.host_removed.fastq.gz"),
    output:
        contigs=join(workpath,"{name}","output","{name}.megahit.contigs.fa"),
        report=join(workpath,"{name}","info","{name}.megahit.contigs_kraken2_report.txt"),
        k2txt=join(workpath,"{name}","kraken2","{name}.megahit.contigs.kraken2"),
        krona=join(workpath,"{name}","kraken2","{name}.megahit.contigs.kraken2.krona"),
        tmp1=join(workpath,"{name}","temp","{name}.megahit_kraken2.txt"),
        ktaxids=join(workpath,"{name}","temp","{name}.megahit_kraken2_taxid.txt"),
        ktmp=join(workpath,"{name}","temp","{name}.megahit_kraken2.temp"),
        knames=join(workpath,"{name}","temp","{name}.megahit_kraken2_names.txt"),
        kcounts_temp=join(workpath,"{name}","temp","{name}.megahit_kraken2_counts.temp"),
        kcounts_txt=join(workpath,"{name}","temp","{name}.megahit_kraken2_counts.txt"),
        kflt=join(workpath,"{name}","temp","{name}.megahit_kraken2_filtered.txt"),
        ktaxlineage=join(workpath,"{name}","temp","{name}.megahit_kraken2_taxlineage.tsv"),
        kviral_unf=join(workpath,"{name}","output","{name}.megahit_kraken2_unfiltered_viraltable.tsv"),
        kviral_flt=join(workpath,"{name}","output","{name}.megahit_kraken2_filtered_viraltable.tsv"),
        kviral_flt_family=join(workpath,"{name}","output","{name}.megahit_kraken2_filtered_viraltable.family-level.tsv"),
        kviral_flt_genus=join(workpath,"{name}","output","{name}.megahit_kraken2_filtered_viraltable.genus-level.tsv"),
        kviral_flt_species=join(workpath,"{name}","output","{name}.megahit_kraken2_filtered_viraltable.species-level.tsv"),
        cat_class=join(workpath,"{name}","CAT","{name}.megahit.contig2classification.txt"),
        cat_names=join(workpath,"{name}","CAT","{name}.megahit.official_names.txt"),
        cat_summary=join(workpath,"{name}","CAT","{name}.megahit.summary.txt"),
        taxids=join(workpath,"{name}","CAT","{name}.megahit.contig_taxids.txt"),
        tmp2=join(workpath,"{name}","temp","{name}.megahit_CAT.txt"),
        tmp3=join(workpath,"{name}","temp","{name}.megahit.kraken2.viral.names.txt"),
        kraken_contigs=join(workpath,"{name}","output","{name}.megahit.kraken2_viral.contigs.fa"),
        tmp4=join(workpath,"{name}","temp","{name}.megahit.CAT.viral.names.txt"),
        cat_contigs=join(workpath,"{name}","output","{name}.megahit.cat_viral.contigs.fa"),
        sam=temp(join(workpath,"{name}","temp","{name}.megahit.sam")),
        bam=temp(join(workpath,"{name}","temp","{name}.megahit.bam")),
        final=join(workpath,"{name}","output","{name}.megahit.bam"),  
    params:
        rname='megahit',
        script=join(workpath, "workflow", "scripts", "collapse_on_taxa.py"),
        filter_length=config['options']['length_filter'],
        taxonkit_db=config['references']['taxonkit_db'],
        viral_db=config['references']['kraken2_viral_db'],
        cat_db=config['references']['CAT_db'],
        cat_tax=config['references']['CAT_taxonomy'],
        cat_dep=config['references']['CAT_diamond'],
        cat_dir=join(workpath,"{name}","CAT"),
        # Building single-end and paired-end options
        # Metahit SE vs PE conditional option: PE='-1', SE='-r'
        megahit_sepe = lambda w: "-1" if nends[w.name] else "-r",
        # Input host removed R2 fastq: PE='-2 input.r2', SE=''
        r2_option = lambda w: "-2 {0}.R2.trim.host_removed.fastq.gz".format(
            join(workpath, w.name, "trim", w.name)
        ) if nends[w.name] else "",
        # Bowtie2 SE vs PE conditional option: PE='-1', SE='-U'
        bowtie2_sepe = lambda w: "-1" if nends[w.name] else "-U",
    threads: int(allocated("threads", "megahit", cluster))
    # conda: config['conda']['CAT']
    # envmodules: 
    #     config['tools']['megahit'], config['tools']['kraken'], config['tools']['kronatools'],
    #     config['tools']['seqtk'], config['tools']['bowtie'], config['tools']['samtools']
    container: config['images']['metavirs']
    shell: """
    set -x
    mkdir -p {wildcards.name}/megahit
    # megahit output directory cannot
    # exist prior to running
    if [ -d "{wildcards.name}/megahit" ]; then
        rm -rf "{wildcards.name}/megahit"
    fi

    # Assemble reads into contigs with megahit
    megahit -t {threads} \\
        {params.megahit_sepe} {input.r1} {params.r2_option} \\
        --out-prefix=megahit \\
        -o {wildcards.name}/megahit
    tr ' ' '_' < {wildcards.name}/megahit/megahit.contigs.fa > {output.contigs}

    # Annotate the contigs with kraken2
    kraken2  --threads {threads} \\
        --db {params.viral_db} {output.contigs} \\
        --report {output.report} \\
    > {output.k2txt}

    awk -v OFS='\\t' '{{split($2,a,"_"); \\
        split(a[4],b,"="); if ($1 == "C") print $2,$3,b[2]}}' \\
        {output.k2txt} > {output.krona}
    cp {output.krona} {output.tmp1}

    # Check the number of assembled contigs
    # prior to running CAT. CAT contigs can
    # fail if this file is empty or only 
    # contains a few contigs.
    mkdir -p {params.cat_dir}/
    n_contigs=$(grep -c '^>' {output.contigs} || true)
    echo "Number of contigs from metaspades: ${{n_contigs}}"
    {{ CAT contigs -n {threads} \\
        --force \\
        -c {output.contigs} \\
        -d {params.cat_db} \\
        -t {params.cat_tax} \\
        --out_prefix {params.cat_dir}/{wildcards.name}.megahit \\
        --path_to_diamond {params.cat_dep}
    CAT add_names -i {output.cat_class} \\
        -o {output.cat_names} \\
        -t {params.cat_tax} \\
        --only_official
    CAT summarise \\
        -c {output.contigs} \\
        -i {output.cat_names} \\
        -o {output.cat_summary}
    }} || {{
        # CAT can fail if provided an empty FASTA file
        echo "WARNING: CAT failed!" 
        echo "This could be due to an issue upstream,"
        echo "such as insufficent sequencing depth,"
        echo "low viral load, or poor quality reads."
        echo "Each of these issues can lead to poor"
        echo "quality assemblies and a low number of"
        echo "assembled contigs."
        echo "Please check the input and output of"
        echo "megahit to troubleshoot the issue!"
        touch {output}
    }}

    # Try to grep for Virsuses, may pipefail
    # if nothing is found, touch output file 
    {{ grep "Viruses:" {output.cat_names} \\
        | awk -v OFS='\\t' '{{split($8,a,";"); \\
            split(a[split($8,a,";")],b,"*"); \\
            split($1,c,"_"); split(c[4],d,"="); \\
            print $1,b[1],d[2]}}' \\
        > {output.taxids} 
    }} || touch {output.taxids}
    cp {output.taxids} {output.tmp2}

    cut -f1 {output.krona} > {output.tmp3}
    seqtk subseq \\
        {output.contigs} \\
        {output.tmp3} \\
    > {output.kraken_contigs}
    cut -f1 {output.taxids} > {output.tmp4}
    seqtk subseq \\
        {output.contigs} \\
        {output.tmp4} \\
    > {output.cat_contigs}

    # Build an index with assembled contigs
    # and align reads against them
    {{ bowtie2-build --threads {threads} \\
        {output.contigs} \\
        {wildcards.name}/temp/{wildcards.name}_megahit
    bowtie2 -p {threads} \\
        -x {wildcards.name}/temp/{wildcards.name}_megahit \\
        {params.bowtie2_sepe} {input.r1} {params.r2_option} \\
        -S {output.sam}
    samtools view -@ {threads} \\
        -bS {output.sam} > {output.bam} 
    samtools sort \\
        -@ {threads} \\
        {output.bam}  \\
        -o {output.final}
    samtools index \\
        -@ {threads} {output.final}
    }} || {{
        # Bowtie failed, maybe due to empty contigs FASTA file
        echo "WARNING: An error occurred while running bowtie!"
        echo "Please see the check the inputs to bowtie2-build"
        echo "and bowtie2 to debug the issue."
        touch {output}
    }} 

    # Filter results based on contig length,
    # get counts for the number of reads aligning
    # to each contig, and create a counts table
    # with counts and taxonomic classifications,
    # Get taxids from intermediate results
    cut -f2 {output.tmp1} \\
    > {output.ktaxids}
    # Get contig names and coverage
    cut -f1,3 {output.tmp1} \\
    > {output.ktmp}
    # Get contig names
    cut -f1 {output.tmp1} \\
    > {output.knames}
    # Get counts of each contig in matching order
    set +o pipefail # grep will fail -f file is empty
    samtools idxstats {output.final} \\
        | grep -f {output.knames} \\
        | cut -f3 \\
    > {output.kcounts_temp}
    set -o pipefail
    # Create text file with contig name, 
    # coverage, and counts
    paste \\
        {output.ktmp} \\
        {output.kcounts_temp} \\
    > {output.kcounts_txt}
    # Create a unfiltered and contig length 
    # filtered table with taxonomic information
    awk -F '\\t' '{{split($1,a,"="); if (a[4] >= {params.filter_length}) print $0}}' \\
        {output.tmp1} \\
    > {output.kflt}
    cat {output.ktaxids} \\
        | taxonkit reformat \\
            --data-dir {params.taxonkit_db} \\
            -I 1 \\
            -r "Unassigned" \\
            -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}\\t{{t}}" \\
    > {output.ktaxlineage}
    # Unfiltered viral taxonomic table
    echo -e "contig\\tcov\\tcount\\ttaxid\\tkingdom\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\tstrain" \\
    > {output.kviral_unf}
    paste \\
        {output.kcounts_txt} \\
        {output.ktaxlineage} \\
    >> {output.kviral_unf}
    # Filtered viral taxonomic table
    head -1 {output.kviral_unf} \\
    > {output.kviral_flt}
    awk -F '\\t' '{{split($1,a,"="); if (NR==1) {{print}}  else if (a[4] >= {params.filter_length}) print $0}}' \\
        {output.kviral_unf} \\
    >> {output.kviral_flt}
    # Collapse and aggregate cov/counts at family-level
    python3 {params.script} \\
        --input {output.kviral_flt} \\
        --output {output.kviral_flt_family} \\
        --collapse-on-taxa family
    # Collapse and aggregate cov/counts at genus-level
    python3 {params.script} \\
        --input {output.kviral_flt} \\
        --output {output.kviral_flt_genus} \\
        --collapse-on-taxa genus
    # Collapse and aggregate cov/counts at species-level
    python3 {params.script} \\
        --input {output.kviral_flt} \\
        --output {output.kviral_flt_species} \\
        --collapse-on-taxa species
    """


rule blast_metaspades_contigs:
    """
    Data-processing step to blast assembled contigs against nt virsuses database.
    For more information see:
    https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html
    https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses-nucl-metadata.json
    https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.0/
    @Input:
        FASTA file of assembled metaspades contigs (scatter)
    @Output:
        Blast text file containing alignment results of metaspades contigs,
    """
    input:
        contigs=join(workpath,"{name}","output","{name}.metaspades.contigs.fa"),
    output:
        blast=join(workpath,"{name}","output","{name}.metaspades_blast.tsv"),
    params:
        rname='blastmetaspades',
        blast_db=config['references']['blast_viral_db'],
        header=join(workpath,"{name}","output","{name}.metaspades_blast.header.tsv"),
        tmp=join(workpath,"{name}","output","{name}.metaspades_blast.tmp.tsv"),
        outfmt_tabs='\\t'.join([
            'qaccver', 'saccver', 'staxid', 
            'ssciname', 'scomname', 'sblastname', 
            'sscinames', 'stitle', 'pident', 'qlen', 
            'length', 'mismatch', 'gapopen', 'qstart',
            'qend', 'sstart', 'send', 'evalue', 'bitscore'
        ]),
        outfmt_spaces=' '.join([
            'qaccver', 'saccver', 'staxid', 
            'ssciname', 'scomname', 'sblastname', 
            'sscinames', 'stitle', 'pident', 'qlen', 
            'length', 'mismatch', 'gapopen', 'qstart',
            'qend', 'sstart', 'send', 'evalue', 'bitscore'
        ]),
    threads: int(allocated("threads", "blast_metaspades_contigs", cluster))
    container: config['images']['blast']
    shell: """
    # BLAST metaspades contigs against viral database
    echo -e "{params.outfmt_tabs}" \\
    > {params.header}
    blastn -query {input.contigs} \\
        -db {params.blast_db} \\
        -out {params.tmp} \\
        -outfmt "6 {params.outfmt_spaces}" \\
        -max_target_seqs 1 \\
        -num_threads {threads}
    # Adding header to blast output
    cat {params.header} \\
        {params.tmp} \\
    > {output.blast}
    rm -f {params.header} {params.tmp}
    """


rule blast_metaspades_xlsx:
    """
    Data-processing step to aggregate the blast-ed metaspades contigs results into
    a single XLSX file. Within the output XLSX file, there will be a tab for each
    sample's blast results. The first tab in the spreadsheet will contain concat
    blast results across all samples.
    @Input:
        Blast texts file containing alignment results of metaspades contigs (gather)
    @Output:
        Blast XLSX file containing alignment results of metaspades contigs for all samples
    """
    input:
        blasts=expand(join(workpath,"{name}","output","{name}.metaspades_blast.tsv"), name=pe_samples),
    output:
        tsv=join(workpath,"Project","metaspades_blast.tsv"), 
        excel=join(workpath,"Project","metaspades_blast.xlsx"),
    params:
        rname='blastmetaspadesxlsx',
        outfmt_tabs='\\t'.join([
            'sample', 'qaccver', 'saccver', 'staxid', 
            'ssciname', 'scomname', 'sblastname', 
            'sscinames', 'stitle', 'pident', 'qlen', 
            'length', 'mismatch', 'gapopen', 'qstart',
            'qend', 'sstart', 'send', 'evalue', 'bitscore'
        ]),
        extension=".metaspades_blast.tsv",
        script=join(workpath, "workflow", "scripts", "files2spreadsheet.py"),
    threads: int(allocated("threads", "blast_metaspades_xlsx", cluster))
    container: config['images']['blast']
    shell: """
    # Create aggregated single file 
    # results of BLAST-ed metaspades 
    # contigs against viral database
    echo -e "{params.outfmt_tabs}" \\
    > {output.tsv}
    # Adding a column to the output
    # that contains each sample name
    awk -F '\\t' -v OFS='\\t' \\
        'FNR>1 {{sub(".*/", "", FILENAME); sub("{params.extension}", "", FILENAME); print FILENAME, $0}}' \\
        {input.blasts} \\
    >> {output.tsv}

    # Create an excel spreadsheet 
    # containing the blast results 
    # of aligning the metaspades
    # contigs against NCBI viral db 
    {params.script} \\
        --add-auto-filters \\
        --rm-suffix "{params.extension}" \\
        --input {output.tsv} {input.blasts} \\
        --output {output.excel}
    """


rule blast_megahit_contigs:
    """
    Data-processing step to blast assembled contigs against nt virsuses database.
    For more information see:
    https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html
    https://ftp.ncbi.nlm.nih.gov/blast/db/nt_viruses-nucl-metadata.json
    https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.0/
    @Input:
        FASTA file of assembled megahit contigs (scatter)
    @Output:
        Blast text file containing alignment results of megahit contigs,
    """
    input:
        contigs=join(workpath,"{name}","output","{name}.megahit.contigs.fa"),
    output:
        blast=join(workpath,"{name}","output","{name}.megahit_blast.tsv"),
    params:
        rname='blastmegahit',
        blast_db=config['references']['blast_viral_db'],
        header=join(workpath,"{name}","output","{name}.megahit_blast.header.tsv"),
        tmp=join(workpath,"{name}","output","{name}.megahit_blast.tmp.tsv"),
        outfmt_tabs='\\t'.join([
            'qaccver', 'saccver', 'staxid', 
            'ssciname', 'scomname', 'sblastname', 
            'sscinames', 'stitle', 'pident', 'qlen', 
            'length', 'mismatch', 'gapopen', 'qstart',
            'qend', 'sstart', 'send', 'evalue', 'bitscore'
        ]),
        outfmt_spaces=' '.join([
            'qaccver', 'saccver', 'staxid', 
            'ssciname', 'scomname', 'sblastname', 
            'sscinames', 'stitle', 'pident', 'qlen', 
            'length', 'mismatch', 'gapopen', 'qstart',
            'qend', 'sstart', 'send', 'evalue', 'bitscore'
        ]),
    threads: int(allocated("threads", "blast_megahit_contigs", cluster))
    container: config['images']['blast']
    shell: """
    # BLAST megahit contigs against viral database
    echo -e "{params.outfmt_tabs}" \\
    > {params.header}
    blastn -query {input.contigs} \\
        -db {params.blast_db} \\
        -out {params.tmp} \\
        -outfmt "6 {params.outfmt_spaces}" \\
        -max_target_seqs 1 \\
        -num_threads {threads}
    # Adding header to blast output
    cat {params.header} \\
        {params.tmp} \\
    > {output.blast}
    rm -f {params.header} {params.tmp}
    """


rule blast_megahit_xlsx:
    """
    Data-processing step to aggregate the blast-ed megahit contigs results into
    a single XLSX file. Within the output XLSX file, there will be a tab for each
    sample's blast results. The first tab in the spreadsheet will contain concat
    blast results across all samples.
    @Input:
        Blast texts file containing alignment results of megahit contigs (gather)
    @Output:
        Blast XLSX file containing alignment results of megahit contigs for all samples
    """
    input:
        blasts=expand(join(workpath,"{name}","output","{name}.megahit_blast.tsv"), name=samples),
    output:
        tsv=join(workpath,"Project","megahit_blast.tsv"), 
        excel=join(workpath,"Project","megahit_blast.xlsx"),
    params:
        rname='blastmegahitxlsx',
        outfmt_tabs='\\t'.join([
            'sample', 'qaccver', 'saccver', 'staxid', 
            'ssciname', 'scomname', 'sblastname', 
            'sscinames', 'stitle', 'pident', 'qlen', 
            'length', 'mismatch', 'gapopen', 'qstart',
            'qend', 'sstart', 'send', 'evalue', 'bitscore'
        ]),
        extension=".megahit_blast.tsv",
        script=join(workpath, "workflow", "scripts", "files2spreadsheet.py"),
    threads: int(allocated("threads", "blast_megahit_xlsx", cluster))
    container: config['images']['blast']
    shell: """
    # Create aggregated single file 
    # results of BLAST-ed megahit 
    # contigs against viral database
    echo -e "{params.outfmt_tabs}" \\
    > {output.tsv}
    # Adding a column to the output
    # that contains each sample name
    awk -F '\\t' -v OFS='\\t' \\
        'FNR>1 {{sub(".*/", "", FILENAME); sub("{params.extension}", "", FILENAME); print FILENAME, $0}}' \\
        {input.blasts} \\
    >> {output.tsv}

    # Create an excel spreadsheet 
    # containing the blast results 
    # of aligning the megahit
    # contigs against NCBI viral db 
    {params.script} \\
        --add-auto-filters \\
        --rm-suffix "{params.extension}" \\
        --input {output.tsv} {input.blasts} \\
        --output {output.excel}
    """


rule krona:
    """
    Reporting step to create a interactive Krona charts based on
    the annotated assembled contigs from the product of of multiple
    assemblers (metaspades and megahit) and annotation tools (CAT and 
    kraken2).
    @Input:
        Annotated contigs from metaspades and megahit (gather per-sample)
    @Output:
        Interactive Krona report for a sample. 
    """
    input:
        get_annotated_contigs
    output:
        report=join(workpath,"{name}","output","{name}.contig.classification.html"),
    params:
        rname='krona',
        krona_ref=config['references']['kronatools'],
    threads: int(allocated("threads", "krona", cluster))
    # envmodules: config['tools']['kraken'], config['tools']['kronatools'],
    container: config['images']['metavirs']
    shell: """
    ktImportTaxonomy \\
        -tax {params.krona_ref} \\
        -m 3 \\
        -o {output.report} \\
        {input}
    """


rule prep_metaquast:
    """
    Data-processing step to prep input for quast. This rule must 
    be logically seperated from metaquast due to the UCSC tool
    dependency. 
    TODO: merge prep_metaquast and metaquast rules together once 
    a docker image has been built with all dependencies.
    @Input:
        Kraken viral aggregrate counts/clade to file (scatter)
    @Output:
        Reformatted text file and split fasta file for metaquast 
    """
    input:
        report=join(workpath,"{name}","info","{name}.reads_kraken2_report.txt"),
    output:
        txt=join(workpath,"{name}","temp","{name}_reads_class_names.txt"),
        tmp=temp(join(workpath,"{name}","temp","{name}.metaquast.tmp.fa")),
        fa=join(workpath,"{name}","temp","{name}.metaquast.fa"),
    params:
        rname='prepmetaq',
        ncbi_viral=config['references']['ncbi_viral_fasta'],
        outdir=join(workpath,"{name}","temp","metaquastref"),
    threads: int(allocated("threads", "prep_metaquast", cluster))
    # envmodules: config['tools']['ucsc']
    container: config['images']['metavirs']
    shell: """
    awk -F '\\t' '{{if ($4 ~ "S") print $6}}' \\
        {input.report} \\
        | sed 's/^ *//g' \\
        | sort \\
        | uniq \\
        | tr ' ' '_' \\
    > {output.txt}

    # Subsets NCBI Viral FASTA file
    # to only include viruses found
    # in the sample in a more fault
    # tolerant manner
    paste - - < {params.ncbi_viral} \\
        > {output.tmp}
    while read default_pattern; do
        # Back pattern to search for
        # if default pattern does not
        # exist in the NCBI viral FASTA
        backup_pattern=$(
            echo "$default_pattern" \\
                | awk -F '_' '{{print $1}}'
        )
        match=$(
            grep "$default_pattern" {output.tmp} \\
                || grep "$backup_pattern" {output.tmp} \\
                || true;
        )
        # Check if there was a match
        # prior to adding to FASTA,
        # avoids adding empty str
        if [ "$match" != "" ]; then  
            echo "$match" \\
                | tr '\\t' '\\n' \\
                | awk -F ',' '{{print $1}}' \\
                | tr '/' '_' \\
                | cut -d '_' -f1-5
        fi
    done < {output.txt} > {output.fa}
    
    mkdir -p {params.outdir}
    faSplit byname {output.fa} {params.outdir}/
    """


rule metaquast:
    """
    Quality-control step of assess the quality of the assembly. 
    @Input:
        Fasta file of assembled contigs from Metaspades and megahit  (scatter)
    @Output:
        Quality-control report of the assembly.
    """
    input:
        contigs=get_contigs_fasta,
        dep=join(workpath,"{name}","temp","{name}.metaquast.fa"),
    output:
        log=join(workpath,"{name}","output","{name}_metaquast","metaquast.log"),
    params:
        rname='metaq',
        ref=join(workpath,"{name}","temp","metaquastref"),
        outdir=join(workpath,"{name}","output","{name}_metaquast"),
    threads: int(allocated("threads", "metaquast", cluster))
    container: config['images']['metavirs']
    shell: """
    metaquast.py \\
        {input.contigs} \\
        -r {params.ref} \\
        --fragmented \\
        --gene-finding \\
        --unique-mapping \\
        -o {params.outdir} \\
        --threads {threads} || {{
        # metaquast can fail if provided an empty FASTA file
        echo "WARNING: metaquast failed!" 
        echo "This could be due to an issue upstream,"
        echo "such as insufficent sequencing depth,"
        echo "low viral load, or poor quality reads."
        echo "Each of these issues can lead to poor"
        echo "quality assemblies and a low number of"
        echo "assembled contigs."
        touch {output}
        exit 0
    }}
    """


rule prep_aggregate_metaspades_scatter:
    """
    Data-processing step to prep input for the aggregated Krona report. Sample
    names are added to text files prior to project-level aggregation. This rule
    is split from prep_aggregate_megahit due to the conditional manner in which 
    metaspades runs for a given sample (does not support single-end data). This
    rule will only run if a sample is paired-end.
    @Input:
        Annotated contigs from metaspades (gather per-sample)
    @Output:
        Annotated contigs with sample names
    """
    input:
        f1=join(workpath,"{name}","temp","{name}.metaspades_kraken2.txt"),
        f2=join(workpath,"{name}","temp","{name}.metaspades_CAT.txt"),
        f3=join(workpath,"{name}","temp","{name}.metaspades_kraken2_filtered.txt"),
    output:
        f1=join(workpath,"Project","temp","{name}.metaspades.contigs.kraken2.krona.meta"),
        f2=join(workpath,"Project","temp","{name}.metaspades.contigs.CAT.krona.meta"),
        f3=join(workpath,"Project","temp","{name}.metaspades_filtered.contigs.kraken2.krona.meta"),
    params:
        rname='prepagg',
        sample='{name}',
    threads: int(allocated("threads", "prep_aggregate_metaspades_scatter", cluster))
    container: config['images']['metavirs']
    shell: """
    # Adding sample name to results per-sample results,
    # metaspades and kraken2
    awk -v var="{params.sample}" \\
        '{{print var"_"$0}}' \\
        {input.f1} \\
    > {output.f1}
    # metaspades and CAT
    awk -v var="{params.sample}" \\
        '{{print var"_"$0}}' \\
        {input.f2} \\
    > {output.f2}
    # Filtered metaspades and kraken
    awk -v var="{params.sample}" \\
        '{{print var"_"$0}}' \\
        {input.f3} \\
    > {output.f3}
    """


rule prep_aggregate_megahit_scatter:
    """
    Data-processing step to prep input for the aggregated Krona report. Sample
    names are added to text files prior to project-level aggregation. This rule
    is split from prep_aggregate_metaspades due to the conditional manner in which 
    megaspades runs for a given sample (does not support single-end data). This
    rule will always run (megahit supports PE and SE data).
    @Input:
        Annotated contigs from megahit (gather per-sample)
    @Output:
        Annotated contigs with sample names
    """
    input:
        f1=join(workpath,"{name}","temp","{name}.megahit_kraken2.txt"),
        f2=join(workpath,"{name}","temp","{name}.megahit_CAT.txt"),
        f3=join(workpath,"{name}","temp","{name}.megahit_kraken2_filtered.txt"),
    output:
        f1=join(workpath,"Project","temp","{name}.megahit.contigs.kraken2.krona.meta"),
        f2=join(workpath,"Project","temp","{name}.megahit.contigs.CAT.krona.meta"),
        f3=join(workpath,"Project","temp","{name}.megahit_filtered.contigs.kraken2.krona.meta"),
    params:
        rname='prepagg',
        sample='{name}',
    threads: int(allocated("threads", "prep_aggregate_megahit_scatter", cluster))
    container: config['images']['metavirs']
    shell: """
    # Adding sample name to results per-sample results,
    # metaspades and kraken2
    awk -v var="{params.sample}" \\
        '{{print var"_"$0}}' \\
        {input.f1} \\
    > {output.f1}
    # metaspades and CAT
    awk -v var="{params.sample}" \\
        '{{print var"_"$0}}' \\
        {input.f2} \\
    > {output.f2}
    # Filtered metaspades and kraken2
    awk -v var="{params.sample}" \\
        '{{print var"_"$0}}' \\
        {input.f3} \\
    > {output.f3}
    """


rule prep_aggregate_metaspades_gather:
    """
    Data-processing step to prep input for the aggregated Krona report.
    All metaspades results are aggregated across all samples. This rule
    is split from pre_aggregate_megahit due to the conditional manner in 
    which metaspades runs for a given sample (does not support single-end 
    data). This rule will only run if a sample is paired-end.
    @Input:
        Annotated contigs from metaspades (gather all-samples)
    @Output:
        Annotated contigs with sample names
    """
    input:
        f1=expand(join(workpath,"Project","temp","{name}.metaspades.contigs.kraken2.krona.meta"), name=pe_samples), 
        f2=expand(join(workpath,"Project","temp","{name}.metaspades.contigs.CAT.krona.meta"), name=pe_samples),
        f3=expand(join(workpath,"Project","temp","{name}.metaspades_filtered.contigs.kraken2.krona.meta"), name=pe_samples),
    output:
        f1=join(workpath,"Project","temp","metaspades_kraken2.krona"),
        f2=join(workpath,"Project","temp","metaspades_CAT.krona"),
        f3=join(workpath,"Project","temp","metaspades_kraken2_filtered.krona"),
    params:
        rname='prepagg',
        search=join(workpath,"Project","temp"),
    threads: int(allocated("threads", "prep_aggregate_metaspades_gather", cluster))
    container: config['images']['metavirs']
    shell: """
    # Create aggregated contig annotations across all samples,
    # Gather metaspades and kraken2 results 
    find {params.search} \\
        -name '*.metaspades.contigs.kraken2.krona.meta' \\
        -exec cat {{}} \\; \\
    > {output.f1}
    # Gather metaspades and CAT results 
    find {params.search} \\
        -name '*.metaspades.contigs.CAT.krona.meta' \\
        -exec cat {{}} \\; \\
    > {output.f2}
    # Gather filtered metaspades and kraken2 results
    find {params.search} \\
        -name '*.metaspades_filtered.contigs.kraken2.krona.meta' \\
        -exec cat {{}} \\; \\
    > {output.f3}
    """


rule prep_aggregate_megahit_gather:
    """
    Data-processing step to prep input for the aggregated Krona report.
    All megahit results are aggregated across all samples. This rule
    is split from pre_aggregate_metaspades due to the conditional manner 
    in which megahit runs for a given sample (does not support single-end 
    data). This rule will only run if a sample is paired-end.
    @Input:
        Annotated contigs from megahit (gather all-samples)
    @Output:
        Annotated contigs with sample names
    """
    input:
        f1=expand(join(workpath,"Project","temp","{name}.megahit.contigs.kraken2.krona.meta"), name=samples), 
        f2=expand(join(workpath,"Project","temp","{name}.megahit.contigs.CAT.krona.meta"), name=samples),
        f3=expand(join(workpath,"Project","temp","{name}.megahit_filtered.contigs.kraken2.krona.meta"), name=samples),
    output:
        f1=join(workpath,"Project","temp","megahit_kraken2.krona"),
        f2=join(workpath,"Project","temp","megahit_CAT.krona"),
        f3=join(workpath,"Project","temp","megahit_kraken2_filtered.krona"),
    params:
        rname='prepagg',
        search=join(workpath,"Project","temp"),
    threads: int(allocated("threads", "prep_aggregate_megahit_gather", cluster))
    container: config['images']['metavirs']
    shell: """
    # Create aggregated contig annotations across all samples,
    # Gather megahit and kraken2 results
    find {params.search} \\
        -name '*.megahit.contigs.kraken2.krona.meta' \\
        -exec cat {{}} \\; \\
    > {output.f1}
    # Gather megahit and CAT results
    find {params.search} \\
        -name '*.megahit.contigs.CAT.krona.meta' \\
        -exec cat {{}} \\; \\
    > {output.f2}
    # Gather filtered megahit and kraken2 results
        find {params.search} \\
        -name '*.megahit_filtered.contigs.kraken2.krona.meta' \\
        -exec cat {{}} \\; \\
    > {output.f3}
    """


rule aggregate:
    """
    Data-processing step to create a project-level, multi-sample aggregated 
    interactive Krona report.
    @Input:
        Annotated contigs with sample names (gather all-samples)
    @Output:
        Multi-sample aggregated Krona report
    """
    input:
        get_aggregated_results,
    output:
        krona=join(workpath,"Project","Project.contig.classification.html"),
    params:
        search=join(workpath,"Project","temp"),
        rname='aggrpt',
        krona_ref=config['references']['kronatools'],
    threads: int(allocated("threads", "aggregate", cluster))
    container: config['images']['metavirs']
    shell: """
    # Generate multi-sample report
    ktImportTaxonomy \\
        -tax {params.krona_ref} \\
        -m 3 \\
        -o {output.krona} \\
        {input}
    """
