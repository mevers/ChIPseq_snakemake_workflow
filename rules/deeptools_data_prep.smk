# Snakemake rules to prepare data for further deeptools processing
#
# Author: Maurits Evers
# License: GPLv3


def get_specific_args(key):
    if key is not None:
        return(key)
    else:
        return("")


# deeptools bamCoverage
rule deeptools_bamCoverage:
    """
    Create bigwig coverage file
    Parameters:
        - Normalise reads per kilobase per million mapped reads
        - Ignore MT X Y rDNA for normalisation
        - Reads are extended to fragment size
        - Duplicate reads are ignored
    """
    version: "1.0"
    conda: "../envs/deeptools.yaml"
    input:
        join(
            config["analysisdir"],
            "{reference_version}/alignment/{sample}.sorted.dedup.bam")
    output:
        join(
            config["analysisdir"],
            "{reference_version}/deeptools/bamCoverage",
            "{sample}.RPKM.bw{binwidth}.bw")
    wildcard_constraints:
        binwidth = "\d+"
    log:
        "logs/deeptools_bamCoverage_{reference_version}_{sample}_RPKM_bw{binwidth}.log"
    params:
        cmd     = "bamCoverage",
        threads = config["deeptools"]["threads"],
        ignore = config["deeptools"]["ignore_for_normalisation"]
    shell:
        """
            {params.cmd} \
            --bam {input} \
            --outFileName {output} \
            --outFileFormat bigwig \
            --binSize {wildcards.binwidth} \
            --normalizeUsingRPKM \
            --ignoreForNormalization {params.ignore} \
            --extendReads \
            -p {params.threads} \
            -v \
            --ignoreDuplicates > {log}
        """


# deeptools bamCompare
rule deeptools_bamCompare:
    """
    Create bigwig file comparing normalised IP vs. control read coverage
    Parameters:
        - Normalise reads per kilobase per million mapped reads
        - Ignore MT X Y rDNA for normalisation
        - Reads are extended to fragment size
        - Duplicate reads are ignored
    """
    version: "1.0"
    conda: "../envs/deeptools.yaml"
    input:
        ChIP = join(
            config["analysisdir"],
            "{reference_version}/alignment",
            "IP_{treatment}_{rep}.sorted.bam"),
        control = join(
            config["analysisdir"],
            "{reference_version}/pooled_control_bam",
            "control_{treatment}_merged.sorted.bam")
    output:
        join(
            config["analysisdir"],
            "{reference_version}/deeptools/bamCompare/{operation}",
            "IP_{treatment}_{rep}_vs_pooled_control.normSES.bw{binwidth}.bw")
    log:
        join(
            "logs",
            "deeptools_bamCompare_{operation}_{reference_version}_" +
            "IP_{treatment}_{rep}_vs_pooled_control.normSES.bw{binwidth}.log")
    threads: 8
    params:
        cmd = "bamCompare",
        threads = config["deeptools"]["threads"],
        ignore = config["deeptools"]["ignore_for_normalisation"],
        other_args = get_specific_args(config["deeptools"]["bamCompare_specific_args"])
    shell:
        """
            {params.cmd} \
            -b1 {input.ChIP} \
            -b2 {input.control} \
            --operation {wildcards.operation} \
            --binSize {wildcards.binwidth} \
            --ignoreForNormalization {params.ignore} \
            {params.other_args} \
            -p {params.threads} \
            -v \
            --ignoreDuplicates \
            -o {output} > {log}
        """


# deeptools computeMatrix for coverage tracks
rule deeptools_computeMatrix_coverage:
    """
    Calculate normalised coverage score per gene
    Parameters:
	    - scale-region mode where all regions are stretched/shrunken to the
          TES->TSS of every gene
	    - Upstream/downstream region defined in config.yaml
        - Skip zero read regions
    """
    version: "1.0"
    conda: "../envs/deeptools.yaml"
    input:
        bw = expand(join(
            config["analysisdir"],
            "{{reference_version}}/deeptools/bamCoverage",
	        "{sample}.RPKM.bw{{binwidth}}.bw"),
            sample = sorted(config["units"].keys())),
        gtf = join(
            config["refdir"],
            config["reference"]["id"],
            config["reference"]["bed_RefSeq"])
    output:
        mat = join(
            config["analysisdir"],
            "{reference_version}/deeptools/computeMatrix",
	        "scoreMatrix.coverage.bw{binwidth}.mat.gz"),
        bed = join(
            config["analysisdir"],
            "{reference_version}/deeptools/computeMatrix",
	        "scoreMatrix.coverage.bw{binwidth}.bed")
    wildcard_constraints:
        binwidth = "\d+"
    log:
        "logs/deeptools_computeMatrix_coverage_{reference_version}_bw{binwidth}.log"
    params:
        cmd        = "computeMatrix",
        upstream   = config["deeptools"]["cm_upstream"],
        downstream = config["deeptools"]["cm_downstream"],
        labels     = expand("{sample}", sample = sorted(config["units"].keys())),
        threads    = config["deeptools"]["threads"],
    shell:
        """
            {params.cmd} scale-regions \
            --scoreFileName {input.bw} \
            --regionsFileName {input.gtf} \
            --outFileName {output.mat} \
            --outFileSortedRegions {output.bed} \
            --upstream {params.upstream} \
            --downstream {params.downstream} \
            --skipZeros \
            --binSize {wildcards.binwidth} \
            --samplesLabel {params.labels} \
            -p {params.threads} \
            --verbose > {log} 2>&1
        """


# deeptools computeMatrix for log2 ratio IP vs. control
rule deeptools_computeMatrix_comparison:
    """
    Calculate comparison score IP vs. control per gene
    Parameters:
	    - scale-region mode where all regions are stretched/shrunken to the
          TES->TSS of every gene
	    - Upstream/downstream region defined in config.yaml
        - Skip zero read regions
    """
    version: "1.0"
    conda: "../envs/deeptools.yaml"
    input:
        bw = [join(
            config["analysisdir"],
            "{reference_version}/deeptools/bamCompare/{operation}",
            re.sub("IP_", "", sample) + ".IP_vs_control.RPKM.bw{binwidth}.bw")
            for sample in sorted(config["units"].keys()) if "IP_" in sample],
        gtf = join(
            config["refdir"],
            config["reference"]["id"],
            config["reference"]["bed_RefSeq"])
    output:
        mat = join(
            config["analysisdir"],
            "{reference_version}/deeptools/computeMatrix",
	        "scoreMatrix.{operation}.bw{binwidth}.mat.gz"),
        bed = join(
            config["analysisdir"],
            "{reference_version}/deeptools/computeMatrix",
	        "scoreMatrix.{operation}.bw{binwidth}.bed")
    wildcard_constraints:
        binwidth = "\d+",
        operation = "(log2|ratio|subtract|add|reciprocal_ratio)"
    log:
        "logs/deeptools_computeMatrix_{operation}_{reference_version}_bw{binwidth}.log"
    params:
        cmd        = "computeMatrix",
        upstream   = config["deeptools"]["cm_upstream"],
        downstream = config["deeptools"]["cm_downstream"],
        labels     = [re.sub("IP_", "", sample)
            for sample in sorted(config["units"].keys()) if "IP_" in sample],
        threads    = config["deeptools"]["threads"],
    shell:
        """
            {params.cmd} scale-regions \
            --scoreFileName {input.bw} \
            --regionsFileName {input.gtf} \
            --outFileName {output.mat} \
            --outFileSortedRegions {output.bed} \
            --upstream {params.upstream} \
            --downstream {params.downstream} \
            --skipZeros \
            --binSize {wildcards.binwidth} \
            --samplesLabel {params.labels} \
            -p {params.threads} \
            --verbose > {log} 2>&1
        """
