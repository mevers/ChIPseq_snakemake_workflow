# Snakemake rules for quality control
#
# Author: Maurits Evers
# License: GPLv3
# Original date: 21-10-2016
# Last changed: 25-07-2017


# Quality control using FASTQC
# Folder/file hierarchy: config["fastqdir"]/{file}
rule process_fastqc_report:
    input:
        join(config["fastqdir"], "{file}.fastq.gz")
    output:
        join(config["fastqcdir"], "{file}_fastqc.zip")
    threads: 1
    params:
        cmd = config["fastqc"]["cmd"],
        out = join(config["fastqcdir"])
    version: "1.0"
    shell:
        """
            {params.cmd} \
            -f fastq \
            -o {params.out} \
            {input}
        """


rule select_polyTA_reads:
    version: "1.0"
    input:
        reads = lambda wildcards: expand(join(
            config["fastqdir"],
            "{file}"),
            file = config["units"][wildcards.sample])
    output:
        r1 = join(
            config["analysisdir"],
            "{reference_version}/processed_data/{sample}_R1.fastq.gz"),
        r2 = join(
            config["analysisdir"],
            "{reference_version}/processed_data/{sample}_R2.fastq.gz")
    params:
        dir = join(config["analysisdir"], "processed_data")
    script:
        "../scripts/select_polyTA_reads.py"
