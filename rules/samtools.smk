# Snakemake rules involving samtools
# http://samtools.sourceforge.net
#
# Note: This rule file should be included in bam.rules.
# Global variables (e.g. directory paths) should be declared in
# bam.rules and are passed down to the individual tool rule files!
#
# Author: Maurits Evers
# License: GPLv3
# Original date: 22-10-2016
# Last changed: 15-08-2017


# Sort and index BAM file
# Note: samtools sort changed its way to specify commandline
# options from version <=0.1.19 to 1.x
# This will potentially break the workflow if run on a machine
# with samtools other than 1.x
rule samtools_sort_and_index:
    input:
        join(
            config["analysisdir"],
            "{reference_version}/alignment/{sample}.bam")
    output:
        bam = join(
            config["analysisdir"],
            "{reference_version}/alignment/{sample}.sorted.bam"),
        bai = join(
            config["analysisdir"],
            "{reference_version}/alignment/{sample}.sorted.bam.bai")
    version: "1.0"
    params:
        cmd = config["samtools"]["cmd"]
    shell:
        """
            {params.cmd} sort -o {output.bam} {input};
            {params.cmd} index {output.bam};
        """


rule samtools_pool_control_and_index:
    input:
        lambda wildcards: expand(join(
            config["analysisdir"],
            "{{reference_version}}/alignment/{sample}.sorted.bam"),
            sample = config["samples"][wildcards.treatment]["control"])
    output:
        join(
            config["analysisdir"],
            "{reference_version}/pooled_control_bam",
            "control_{treatment}_merged.sorted.bam")
    version: "1.0"
    params:
        cmd = config["samtools"]["cmd"]
    shell:
        """
            {params.cmd} merge {output} {input};
            {params.cmd} index {output};
        """


# Remove duplicates using samtools rmdup
rule samtools_rmdup_and_index:
    input:
        join(
            config["analysisdir"],
            "{reference_version}/picard-tools/dupes/{sample}_markedDupes.bam")
    output:
        bam = join(
            config["analysisdir"],
            "{reference_version}/alignment/{sample}.sorted.dedup.bam"),
        bai = join(
            config["analysisdir"],
            "{reference_version}/alignment/{sample}.sorted.dedup.bam.bai")
    params:
        cmd = config["samtools"]["cmd"]
    version: "1.0"
    shell:
        """
            {params.cmd} rmdup {input} {output.bam};
            {params.cmd} index {output.bam}
        """


# Calculate and store statistics using samtools flagstat
rule flagstat_bam:
    input:
        join(
            config["analysisdir"],
            "{reference_version}/alignment/{sample}.bam")
    output:
        join(
            config["analysisdir"],
            "{reference_version}/samtools/flagstat/flagstat_{sample}.txt")
    params:
        cmd = config["samtools"]["cmd"]
    version: "1.0"
    shell:
        """
            {params.cmd} flagstat {input} > {output}
        """


# Extract duplicate reads flagged by picard-tools MarkDuplicates
rule extract_read_dupes:
    input:
        join(config["analysisdir"], config["reference"]["id"],
            "picard-tools/dupes/{sample}_markedDupes.bam")
    output:
        join(config["analysisdir"], config["reference"]["id"],
            "picard-tools/dupes/{sample}.dupes.bam"),
    params:
        cmd = config["samtools"]["cmd"]
    version: "1.0"
    shell:
        """
            {params.cmd} view -bf 0x400 {input} > {output}
        """
