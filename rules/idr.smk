# Snakemake rules for irreproducible discovery rate (IDR) analysis
#
# Author: Maurits Evers
# License: GPLv3


# idr
rule idr:
    """
    Perform irreproducibile discovery analysis
    Tested on IDR 2.0.3 from https://github.com/nboley/idr
    """
    version: "1.0"
    input:
        peak1 = join(
            config["analysisdir"],
            "{reference_version}/macs2/callpeak",
            "IP_{treatment}_{rep1}_vs_pooled_control_peaks.narrowPeak"),
        peak2 = join(
            config["analysisdir"],
            "{reference_version}/macs2/callpeak",
            "IP_{treatment}_{rep2}_vs_pooled_control_peaks.narrowPeak")
    output:
        join(
            config["analysisdir"],
            "{reference_version}/idr/",
            "{treatment}_IDR_{rep1}_vs_{rep2}.narrowPeak")
    log:
        "logs/idr_{reference_version}_{treatment}_{rep1}_vs_{rep2}.log"
    params:
        cmd = "idr"
    shell:
        """
            {params.cmd} \
            --samples {input} \
            --input-file-type narrowPeak \
            --output-file {output} \
            --plot \
            --peak-merge-method avg \
            --verbose > {log} 2>&1
        """


rule idr_plus:
    """
    Perform irreproducibile discovery analysis for >2 samples
    Tested on IDR 2.0.3 from https://github.com/nboley/idr
    """
    version: "1.0"
    input:
        peaks = lambda wildcards: expand(join(
            config["analysisdir"],
            "{{reference_version}}/macs2/callpeak",
            "{sample}_vs_pooled_control_peaks.narrowPeak"),
            sample = config["samples"][wildcards.treatment]["IP"])
    output:
        join(
            config["analysisdir"],
            "{reference_version}/idr/",
            "{treatment}_IDR_summarised.narrowPeak")
    log:
        "logs/idr_{reference_version}_{treatment}.log"
    script:
        "../scripts/idr_plus.py"
