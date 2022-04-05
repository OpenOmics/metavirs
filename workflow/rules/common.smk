# Rules common to pipeline, irrespective if the data is single-end or paired-end
from scripts.common import abstract_location, allocated


rule fc_lane:
    """
    Quality-control step to get flowcell and lane information from FastQ file.
    FastQ files generated with older versions of Casava or downloaded from
    SRA have a different format than newer FastQ files generated with the
    current version of Casava. It is worth noting that FastQ files downloaded from SRA
    or FastQ files generated with Casava version < 1.8 do not have Flowcell
    IDs in its sequence indentifer. If a FastQ file does not have Flowcell IDs,
    the Machine or Instrument ID is grabbed instead.
    @Input:
        Raw FastQ R1 file (scatter)
    @Output:
        Text file containing information about the FastQ file
    """
    input:
        r1=join(workpath,"{name}.R1.fastq.gz"),
    output:
        info=join(workpath,"{name}","rawQC","{name}.fastq.info.txt")
    params:
        rname='fc_lane',
        get_flowcell_lanes=join("workflow", "scripts", "get_flowcell_lanes.py"),
    threads: int(allocated("threads", "fc_lane", cluster))
    # envmodules: config['tools']['python']
    container: config['images']['metavirs']
    shell: """
    python {params.get_flowcell_lanes} \\
        {input.r1} \\
        {wildcards.name} > {output.info}
    """
