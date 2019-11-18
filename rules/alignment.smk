# Snakemake rules for read alignment
#
# Author: Maurits Evers
# License: GPLv3
# Original date: 21-10-2016
# Last changed: 16-08-2017


# Read alignment using bowtie2
# Note: We discard unmapped reads.
# Note: We automatically decide on whether we have SE or PE reads
# based on the length of input.reads
rule bowtie2_align_reads:
    input:
        idx = IDX,
        reads = lambda wildcards: expand(join(
            config["fastqdir"],
            "{file}"),
            file = config["units"][wildcards.sample]),
        qc = lambda wildcards: expand(join(
            config["fastqcdir"],
            "{file}"),
            file = [w.replace(".fastq.gz", "_fastqc.zip")
                for w in config["units"][wildcards.sample]])
    output:
        join(
            config["analysisdir"],
            "{reference_version}/alignment/{sample}.bam")
    log:
        "logs/bowtie2_{reference_version}_{sample}.log"
    params:
        cmd      = config["bowtie2"]["cmd"],
        ref      = re.sub(".fa", "", REF),
        in_fmt   = config["bowtie2"]["in_fmt"],
        phred    = config["bowtie2"]["phred"],
        maxins   = config["bowtie2"]["maxins"],
        mismatch = config["bowtie2"]["mismatch"],
        threads  = config["bowtie2"]["threads"]
    version: "1.0"
    run:
        if (len(input.reads) == 2):
            shell(" \
                {params.cmd} \
                {params.phred} \
                --no-mixed \
                --no-discordant \
                --maxins {params.maxins} \
                -N {params.mismatch} \
                --threads {params.threads} \
                -x {params.ref} \
                -1 {input.reads[0]} \
                -2 {input.reads[1]} \
                2> {log} \
                | samtools view -bS -F4 - > {output} \
            ")
        elif (len(input.reads) == 1):
            shell(" \
                {params.cmd} \
                {params.phred} \
                -N {params.mismatch} \
                --threads {params.threads} \
                -x {params.ref} \
                -U {input.reads} \
                2> {log} \
                | samtools view -bS -F4 - > {output} \
            ")
