# Snakemake rules to generate deeptools plots
#
# Author: Maurits Evers
# License: GPLv3


# deeptools computeMatrix for log2 ratio IP vs. control
rule deeptools_metagene_heatmap:
    """
    Plot metagene heatmap of scores generated from computeMatrix
    Parameters:
	    - scale-region mode where all regions are stretched/shrunken to the
          TES->TSS of every gene
	    - Upstream/downstream region defined in config.yaml
        - Skip zero read regions
    """
    version: "1.0"
    conda: "../envs/deeptools.yaml"
    input:
        mat = join(
            config["analysisdir"],
            "{reference_version}/deeptools/computeMatrix",
            "scoreMatrix.{operation}.bw{binwidth}.mat.gz")
    output:
        pdf = join(
            config["analysisdir"],
            "{reference_version}/deeptools/plotHeatmap",
            "metageneplot.{operation}.bw{binwidth}.pdf"),
        bed = join(
            config["analysisdir"],
            "{reference_version}/deeptools/plotHeatmap",
            "metageneplot.{operation}.bw{binwidth}.bed")
    wildcard_constraints:
        binwidth = "\d+"
    log:
        "logs/deeptools_plotHeatmap_{reference_version}_{operation}_bw{binwidth}.log"
    params:
        cmd        = "plotHeatmap",
        threads    = config["deeptools"]["threads"],
    shell:
        """
            {params.cmd} \
            --matrixFile {input.mat} \
            --outFileName {output.pdf} \
            --outFileSortedRegions {output.bed} \
            --verbose > {log} 2>&1
        """
