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
        log1=join(workpath,"{name}","rawQC","{name}.validated.R1.fastq.log"),
        log2=join(workpath,"{name}","rawQC","{name}.validated.R2.fastq.log"),
    params:
        rname='validfq',
        outdir=join(workpath,"{name}","rawQC"),
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
        Raw FastQ file (scatter)
    @Output:
        FastQC report and zip file containing data quality information
    """
    input:
        join(workpath,"{name}.R1.fastq.gz"),
        join(workpath,"{name}.R2.fastq.gz"),
    output:
        join(workpath,"{name}","rawQC","{name}.R1_fastqc.zip"),
        join(workpath,"{name}","rawQC","{name}.R2_fastqc.zip"),
    params:
        rname='rawfqc',
        outdir=join(workpath,"{name}","rawQC"),
    threads: int(allocated("threads", "rawfastqc", cluster))
    # envmodules: config['tools']['fastqc']
    container: config['images']['metavirs']
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
        r1=join(workpath,"{name}","trim","{name}.R1.trim.fastq"),
        r2=join(workpath,"{name}","trim","{name}.R2.trim.fastq"),
    params:
        rname='trimfq',
        adapters=config['references']['adapters'],
    threads: int(allocated("threads", "trim", cluster))
    # envmodules: config['tools']['cutadapt']
    container: config['images']['metavirs']
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
        Trimmed FastQ file (scatter)
    @Output:
        FastQC report and zip file containing data quality information
    """
    input:
        join(workpath,"{name}","trim","{name}.R1.trim.fastq"),
        join(workpath,"{name}","trim","{name}.R2.trim.fastq"),
    output:
        join(workpath,"{name}","QC","{name}.R1.trim_fastqc.zip"),
        join(workpath,"{name}","QC","{name}.R2.trim_fastqc.zip"),
    params:
        rname='fqc',
        outdir=join(workpath,"{name}","QC"),
    threads: int(allocated("threads", "fastqc", cluster))
    # envmodules: config['tools']['fastqc']
    container: config['images']['metavirs']
    shell: """
    fastqc {input} \\
        -t {threads} \\
        -o {params.outdir}
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
        r2=join(workpath,"{name}","trim","{name}.R2.trim.fastq"),
    output:
        sam=temp(join(workpath,"{name}","temp","{name}.host_contaminated.sam")),
        bam=temp(join(workpath,"{name}","temp","{name}.host_contaminated.bam")),
        sorted_bam=temp(join(workpath,"{name}","temp","{name}.host_contaminated.sorted.bam")),
        r1=join(workpath,"{name}","trim","{name}.R1.trim.host_removed.fastq.gz"),
        r2=join(workpath,"{name}","trim","{name}.R2.trim.host_removed.fastq.gz"),
    params:
        rname='rmhost',
        host_index=join(config['references']['host_bowtie2_index'], 'Hosts')
    threads: int(allocated("threads", "remove_host", cluster))
    # envmodules: config['tools']['bowtie'], config['tools']['samtools']
    container: config['images']['metavirs']
    shell: """
    # Align against host to remove contamination
    bowtie2 -p {threads} \\
        -x {params.host_index} \\
        -1 {input.r1} \\
        -2 {input.r2} \\
        -S {output.sam}
    samtools view -@ {threads} \\
        -bS -f 12 -F 256 \\
        {output.sam} > {output.bam}
    samtools sort -n -@ {threads} \\
        {output.bam} \\
        -o {output.sorted_bam}
    samtools fastq -@ {threads} \\
        {output.sorted_bam} \\
        -1 {output.r1} \\
        -2 {output.r2} \\
        -0 /dev/null -s /dev/null -n
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
        r2=join(workpath,"{name}","trim","{name}.R2.trim.host_removed.fastq.gz"),
    output:
        report=join(workpath,"{name}","info","{name}.reads_kraken2_report.txt"),
        k2txt=join(workpath,"{name}","kraken2","{name}.reads.kraken2"),
        kronatxt=join(workpath,"{name}","kraken2","{name}.reads.kraken2.krona"),
        html=join(workpath,"{name}","kraken2","{name}.reads.krona.html"),
    params:
        rname='krakenviral',
        viral_db=config['references']['kraken2_viral_db'],
        krona_ref=config['references']['kronatools'],
    threads: int(allocated("threads", "kraken_viral", cluster))
    # envmodules: config['tools']['kraken'], config['tools']['kronatools']
    container: config['images']['metavirs']
    shell: """
    # Run kraken against viral database
    kraken2 --threads {threads} \\
        --db {params.viral_db} \\
        --paired {input.r1} {input.r2} \\
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
        r2=join(workpath,"{name}","trim","{name}.R2.trim.host_removed.fastq.gz"),
    output:
        contigs=join(workpath,"{name}","output","{name}.metaspades.contigs.fa"),
        report=join(workpath,"{name}","info","{name}.metaspades.contigs_kraken2_report.txt"),
        k2txt=join(workpath,"{name}","kraken2","{name}.metaspades.contigs.kraken2"),
        krona=join(workpath,"{name}","kraken2","{name}.metaspades.contigs.kraken2.krona"),
        tmp1=join(workpath,"{name}","temp","{name}.metaspades_kraken2.txt"),
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
        viral_db=config['references']['kraken2_viral_db'],
        cat_db=config['references']['CAT_db'],
        cat_tax=config['references']['CAT_taxonomy'],
        cat_dep=config['references']['CAT_diamond'],
        cat_dir=join(workpath,"{name}","CAT"),
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
    metaspades.py -t {threads} \\
        -m 240 \\
        -1 {input.r1} \\
        -2 {input.r2} \\
        -o {wildcards.name}/metaspades
    mkdir -p {wildcards.name}/output
    mv {wildcards.name}/metaspades/contigs.fasta {output.contigs}

    kraken2  --threads {threads} \\
        --db {params.viral_db} {output.contigs} \\
        --report {output.report} \\
    > {output.k2txt}
    awk -v OFS='\\t' '{{split($2,a,"_"); if ($1 == "C") print $2,$3,a[6]}}' \\
        {output.k2txt} > {output.krona}
    cp {output.krona} {output.tmp1}

    mkdir -p {params.cat_dir}
    CAT contigs -n {threads} \\
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

    bowtie2-build --threads {threads} \\
        {output.contigs} \\
        {wildcards.name}/temp/{wildcards.name}_metaspades
    bowtie2 -p {threads} \\
        -x {wildcards.name}/temp/{wildcards.name}_metaspades \\
        -1 {input.r1} \\
        -2 {input.r2} \\
        -S {output.sam}
    
    samtools view -@ {threads} \\
        -bS {output.sam} > {output.bam} 
    samtools sort \\
        -@ {threads} \\
        {output.bam}  \\
        -o {output.final}
    samtools index \\
        -@ {threads} {output.final}
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
        r2=join(workpath,"{name}","trim","{name}.R2.trim.host_removed.fastq.gz"),
    output:
        contigs=join(workpath,"{name}","output","{name}.megahit.contigs.fa"),
        report=join(workpath,"{name}","info","{name}.megahit.contigs_kraken2_report.txt"),
        k2txt=join(workpath,"{name}","kraken2","{name}.megahit.contigs.kraken2"),
        krona=join(workpath,"{name}","kraken2","{name}.megahit.contigs.kraken2.krona"),
        tmp1=join(workpath,"{name}","temp","{name}.megahit_kraken2.txt"),
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
        viral_db=config['references']['kraken2_viral_db'],
        cat_db=config['references']['CAT_db'],
        cat_tax=config['references']['CAT_taxonomy'],
        cat_dep=config['references']['CAT_diamond'],
        cat_dir=join(workpath,"{name}","CAT")
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
    megahit -t {threads} \\
        -1 {input.r1} \\
        -2 {input.r2} \\
        --out-prefix=megahit \\
        -o {wildcards.name}/megahit
    tr ' ' '_' < {wildcards.name}/megahit/megahit.contigs.fa > {output.contigs}

    kraken2  --threads {threads} \\
        --db {params.viral_db} {output.contigs} \\
        --report {output.report} \\
    > {output.k2txt}

    awk -v OFS='\\t' '{{split($2,a,"_"); \\
        split(a[4],b,"="); if ($1 == "C") print $2,$3,b[2]}}' \\
        {output.k2txt} > {output.krona}
    cp {output.krona} {output.tmp1}

    mkdir -p {params.cat_dir}/
    CAT contigs -n {threads} \\
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

    bowtie2-build --threads {threads} \\
        {output.contigs} \\
        {wildcards.name}/temp/{wildcards.name}_megahit
    bowtie2 -p {threads} \\
        -x {wildcards.name}/temp/{wildcards.name}_megahit \\
        -1 {input.r1} \\
        -2 {input.r2} \\
        -S {output.sam}
    
    samtools view -@ {threads} \\
        -bS {output.sam} > {output.bam} 
    samtools sort \\
        -@ {threads} \\
        {output.bam}  \\
        -o {output.final}
    samtools index \\
        -@ {threads} {output.final}
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
        f1=join(workpath,"{name}","temp","{name}.metaspades_CAT.txt"),
        f2=join(workpath,"{name}","temp","{name}.metaspades_kraken2.txt"),
        f3=join(workpath,"{name}","temp","{name}.megahit_CAT.txt"),
        f4=join(workpath,"{name}","temp","{name}.megahit_kraken2.txt"),
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
        {input.f1} \\
        {input.f2} \\
        {input.f3} \\
        {input.f4}
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
    
    paste - - < {params.ncbi_viral} \\
        | grep -f {output.txt} \\
        | tr '\\t' '\\n' \\
        | awk -F ',' '{{print $1}}' \\
        | tr '/' '_' \\
        | cut -d '_' -f1-5 \\
    > {output.fa}
    
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
        metaspades=join(workpath,"{name}","output","{name}.metaspades.contigs.fa"),
        megahit=join(workpath,"{name}","output","{name}.megahit.contigs.fa"),
        dep=join(workpath,"{name}","temp","{name}.metaquast.fa"),
    output:
        report=join(workpath,"{name}","output","{name}_metaquast","report.html")
    params:
        rname='metaq',
        ref=join(workpath,"{name}","temp","metaquastref"),
        outdir=join(workpath,"{name}","output","{name}_metaquast"),
    threads: int(allocated("threads", "metaquast", cluster))
    container: config['images']['metavirs']
    shell: """
    metaquast.py \\
        {input.metaspades} \\
        {input.megahit} \\
        -r {params.ref} \\
        --fragmented \\
        --gene-finding \\
        --unique-mapping \\
        -o {params.outdir} \\
        --threads {threads}
    """


rule prep_aggregate:
    """
    Data-processing step to prep input for the aggregated Krona report. Sample
    names are added to text files prior to project-level aggregation.
    @Input:
        Annotated contigs from metaspades and megahit (gather per-sample)
    @Output:
        Annotated contigs with sample names
    """
    input:
        f1=join(workpath,"{name}","temp","{name}.metaspades_kraken2.txt"),
        f2=join(workpath,"{name}","temp","{name}.megahit_kraken2.txt"),
        f3=join(workpath,"{name}","temp","{name}.metaspades_CAT.txt"),
        f4=join(workpath,"{name}","temp","{name}.megahit_CAT.txt"),
    output:
        f1=join(workpath,"Project","temp","{name}.metaspades.contigs.kraken2.krona.meta"),
        f2=join(workpath,"Project","temp","{name}.megahit.contigs.kraken2.krona.meta"),
        f3=join(workpath,"Project","temp","{name}.metaspades.contigs.CAT.krona.meta"),
        f4=join(workpath,"Project","temp","{name}.megahit.contigs.CAT.krona.meta"),
    params:
        rname='prepagg',
        sample='{name}',
    threads: int(allocated("threads", "prep_aggregate", cluster))
    container: config['images']['metavirs']
    shell: """
    awk -v var="{params.sample}" \\
        '{{print var"_"$0}}' \\
        {input.f1} \\
    > {output.f1}
    awk -v var="{params.sample}" \\
        '{{print var"_"$0}}' \\
        {input.f2} \\
    > {output.f2}
    awk -v var="{params.sample}" \\
        '{{print var"_"$0}}' \\
        {input.f3} \\
    > {output.f3}
    awk -v var="{params.sample}" \\
        '{{print var"_"$0}}' \\
        {input.f4} \\
    > {output.f4}
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
        f1=expand(join(workpath,"Project","temp","{name}.metaspades.contigs.kraken2.krona.meta"), name=samples), 
        f2=expand(join(workpath,"Project","temp","{name}.megahit.contigs.kraken2.krona.meta"), name=samples),
        f3=expand(join(workpath,"Project","temp","{name}.metaspades.contigs.CAT.krona.meta"), name=samples),
        f4=expand(join(workpath,"Project","temp","{name}.megahit.contigs.CAT.krona.meta"), name=samples),
    output:
        f1=join(workpath,"Project","temp","metaspades_kraken2.krona"),
        f2=join(workpath,"Project","temp","megahit_kraken2.krona"),
        f3=join(workpath,"Project","temp","metaspades_CAT.krona"),
        f4=join(workpath,"Project","temp","megahit_CAT.krona"),
        krona=join(workpath,"Project","Project.contig.classification.html"),
    params:
        search=join(workpath,"Project","temp"),
        rname='aggrpt',
        krona_ref=config['references']['kronatools'],
    threads: int(allocated("threads", "aggregate", cluster))
    container: config['images']['metavirs']
    shell: """
    # Create aggregated contig annotations
    find {params.search} \\
        -name '*.metaspades.contigs.kraken2.krona.meta' \\
        -exec cat {{}} \\; \\
    > {output.f1}
    find {params.search} \\
        -name '*.megahit.contigs.kraken2.krona.meta' \\
        -exec cat {{}} \\; \\
    > {output.f2}
    find {params.search} \\
        -name '*.metaspades.contigs.CAT.krona.meta' \\
        -exec cat {{}} \\; \\
    > {output.f3}
    find {params.search} \\
        -name '*.megahit.contigs.CAT.krona.meta' \\
        -exec cat {{}} \\; \\
    > {output.f4}
    # Generate multi-sample report
    ktImportTaxonomy \\
        -tax {params.krona_ref} \\
        -m 3 \\
        -o {output.krona} \\
        {output.f1} \\
        {output.f2} \\
        {output.f3} \\ 
        {output.f4}
    """