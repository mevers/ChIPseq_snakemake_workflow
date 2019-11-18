# Snakemake rules for ChIP-seq peak calling
#
# Author: Maurits Evers
# License: GPLv3
# Original date: 23-10-2016
# Last changed: 16-10-2019


rule macs2_predict_fragment_size:
    """
    MACS2 predict fragment size for single-end data
    """
    version: "1.0"
    conda: "../envs/macs2.yaml"
    input:
        ChIP = join(
            config["analysisdir"],
            "{reference_version}/alignment/IP_{treatment}_{rep}.sorted.bam"),
    output:
        join(
            config["analysisdir"],
            "{reference_version}/macs2/predictd",
            "IP_{treatment}_{rep}.R")
    log:
        join(
            "logs",
            "macs2_predictd_{reference_version}_{treatment}_{rep}.log")
    params:
        cmd = "macs2",
        format = config["macs2"]["format"],
        gsize = config["macs2"]["gsize"],
        qvalue = config["macs2"]["qvalue"],
        name = "IP_{treatment}_{rep}.R",
        outdir = join(
            config["analysisdir"],
            "{reference_version}/macs2/predictd")
    shell:
        """
            {params.cmd} predictd \
            -i {input.ChIP} \
            -f {params.format} \
            -g {params.gsize} \
            -m 5 500 \
            --outdir {params.outdir} \
            --rfile {params.name} &> {log} && \
            Rscript -e 'setwd("{params.outdir}"); source("{params.name}")'
        """


rule macs2_callpeak_with_pooled_control:
    """
    MACS2 peak calling
    Parameters:
        - Remove duplicate reads
        - q value cutoff from config["macs2"]["qvalue"]
        - Select regions with enrichment scores in range [5, 500]
        - Reanalyse shape profile to deconvolve subpeaks within each peak
        - format from config["macs2"]["format"]
    """
    version: "1.0"
    conda: "../envs/macs2.yaml"
    input:
        ChIP = join(
            config["analysisdir"],
            "{reference_version}/alignment/IP_{treatment}_{rep}.sorted.bam"),
        ctrl = lambda wildcards: expand(join(
            config["analysisdir"],
            "{{reference_version}}/alignment/{input_sample}.sorted.bam"),
            input_sample = config["samples"][wildcards.treatment]["control"]),
        fragment_d = join(
            config["analysisdir"],
            "{reference_version}/macs2/predictd",
            "IP_{treatment}_{rep}.R")
    output:
        join(
            config["analysisdir"],
            "{reference_version}/macs2/callpeak",
            "IP_{treatment}_{rep}_vs_pooled_control_peaks.narrowPeak")
    log:
        join(
            "logs",
            "macs2_callpeak_{reference_version}_{treatment}_{rep}.log")
    params:
        cmd = "macs2",
        format = config["macs2"]["format"],
        gsize = config["macs2"]["gsize"],
        qvalue = config["macs2"]["qvalue"],
        name = "IP_{treatment}_{rep}_vs_pooled_control",
        outdir = join(
            config["analysisdir"],
            "{reference_version}/macs2/callpeak")
    shell:
        """
            {params.cmd} callpeak \
            -t {input.ChIP} \
            -c {input.ctrl} \
            -f {params.format} \
            -g {params.gsize} \
            --outdir {params.outdir} \
            -n {params.name} \
            -m 5 500 \
            --call-summits \
            -q {params.qvalue} &> {log} \
        """


rule filter_narrowPeak:
    """
    Filter MACS2-based narrowPeak files
    Parameters:
        - Fold enrichment threshold [default: 10]
    """
    version: "1.0"
    input:
        join(
            config["analysisdir"],
            "{reference_version}/macs2/callpeak",
            "{name}.narrowPeak")
    output:
        join(
            config["analysisdir"],
            "{reference_version}/macs2/callpeak",
            "EFfiltered/{name}.narrowPeak")
    log:
        join(
            "logs",
            "filter_narrowPeak_{reference_version}_{name}.log")
    params:
        thresholdEF = 5
    script:
        "../scripts/filter_narrowPeak.py"
