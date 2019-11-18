# Snakemake rules for irreproducible discovery rate (IDR) analysis
#
# Author: Maurits Evers
# License: GPLv3


# idr
rule idr:
    """
    Perform irreproducibile discovery analysis
    """
    version: "1.0"
    input:
        lambda wildcards: [join(
            config["analysisdir"],
            "{reference_version}/macs2/callpeak",
            name + "_vs_pooled_control_peaks.narrowPeak")
            for name in config["samples"][wildcards.treatment]["IP"]]
    output:
        join(
            config["analysisdir"],
            "{reference_version}/idr/",
            "{treatment}_IDRsummarised.narrowPeak")
    log:
        "logs/idr_{reference_version}_{treatment}.log"
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


#
