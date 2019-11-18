# Snakemake rules to process BAM files
#
# Author: Maurits Evers
# License: GPLv3
# Original date: 22-10-2016
# Last changed: 15-08-2017


# Include individual tool rule files
include: "samtools.smk"
include: "bedtools.rules"
include: "picard-tools.rules"
include: "qualimap.rules"


# deeptools multiBamSummary
rule deeptools_multiBamSummary_genome:
    input:
        expand(join(
            config["analysisdir"],
            "{{reference_version}}/alignment/{sample}.sorted.dedup.bam"),
            sample = config["units"].keys())
    output:
        npz = join(
            config["analysisdir"],
            "{reference_version}/deeptools/cov.all.genome.bw{binwidth}.npz"),
        tab = join(
            config["analysisdir"],
            "{reference_version}/deeptools/counts.all.genome.bw{binwidth}.tab")
    params:
        cmd     = "multiBamSummary",
        labels  = expand("{sample}", sample = config["units"].keys()),
        threads = config["deeptools"]["threads"]
    version: "1.0"
    shell:
        """
            {params.cmd} bins \
            --bamfiles {input} \
            -out {output.npz} \
            --outRawCounts {output.tab} \
            --labels {params.labels} \
            --binSize {wildcards.binwidth} \
            -p {params.threads}
        """


# deeptools multiBamSummary peaks only
rule deeptools_multiBamSummary_peaks:
    input:
        BAM = expand(join(
            config["analysisdir"],
            "{{reference_version}}/alignment/{sample}.sorted.dedup.bam"),
            sample = config["units"].keys()),
        BED = join(
            config["analysisdir"],
            config["reference"]["id"],
            "{reference_version}/macs2/merged_peaks.bed")
    output:
        join(
            config["analysisdir"],
            "{reference_version}/deeptools/multiBamSummary_results_peaks.npz")
    params:
        cmd    = "multiBamSummary",
        labels = expand("{sample}", sample = config["units"].keys())
    version: "1.0"
    shell:
        """
            {params.cmd} BED-file \
            --BED {input.BED} \
            --bamfiles {input.BAM} \
            -out {output} \
            --labels {params.labels}
        """


# deeptools plotCorrelation
rule deeptools_plotCorrelation:
    input:
        join(
            config["analysisdir"],
            "{reference_version}/deeptools/cov.all.genome.bw{binwidth}.npz")
    output:
        join(
            config["analysisdir"],
            "{reference_version}/deeptools/heatmap.cor.genome.bw{binwidth}.pdf")
    wildcard_constraints:
        binwidth = "\d+"
    params:
        cmd = "plotCorrelation"
    version: "1.0"
    shell:
        """
            {params.cmd} \
            --corData {input} \
            --plotFile {output} \
            --corMethod spearman \
            --whatToPlot heatmap \
            --skipZeros \
            --plotTitle "Spearman correlation of read counts across genome (binwidth {wildcards.binwidth})" \
            --removeOutliers \
            --plotNumbers \
            --colorMap RdBu \
            --zMin -1 \
            --zMax +1
        """


# deeptools plotPCA
rule deeptools_plotPCA:
    input:
        join(
            config["analysisdir"],
            "{reference_version}/deeptools/cov.all.genome.bw{binwidth}.npz")
    output:
        join(
            config["analysisdir"],
            "{reference_version}/deeptools/PCA.genome.bw{binwidth}.pdf")
    wildcard_constraints:
        binwidth = "\d+"
    params:
        cmd = "plotPCA"
    version: "1.0"
    shell:
        """
            {params.cmd} \
            --corData {input} \
            --plotFile {output} \
            --rowCenter \
            --plotTitle "PCA based on read counts across genome (binwidth {wildcards.binwidth})"
        """


# deeptools plotFingerprint
rule deeptools_plotFingerprint:
    input:
        expand(join(
            config["analysisdir"],
            "{{reference_version}}/alignment{sample}.sorted.dedup.bam"),
            sample = config["units"].keys())
    output:
        pdf = join(
            config["analysisdir"],
            "{reference_version}/deeptools/fingerprint.bw{binwidth}.skipZeros.pdf"),
        tab = join(
            config["analysisdir"],
            "{reference_version}/deeptools/fingerprint.bw{binwidth}.skipZeros.tab")
    wildcard_constraints:
        binwidth = "\d+"
    params:
        cmd     = "plotFingerprint",
        labels  = expand("{sample}", sample = config["units"].keys()),
        threads = config["deeptools"]["threads"]
    version: "1.0"
    shell:
        """
            {params.cmd} \
            --bamfiles {input} \
            --plotFile {output.pdf} \
            --ignoreDuplicates \
            --labels {params.labels} \
            --binSize {wildcards.binwidth} \
            --skipZeros \
            --plotTitle "Fingerprints of different samples" \
            --outRawCounts {output.tab} \
            -p {params.threads}
        """
