# Unused explicit snakemake targets for rule all

# BAM targets: sorted and sorted+deduped BAM files plus indices
ALL_BAM = expand(join(
    config["analysisdir"],
    config["reference"]["id"],
    "alignment/{sample}.{suf}"),
    sample = config["units"].keys(),
    suf = [
        "bam",
        "sorted.bam",
        "sorted.bam.bai",
        "sorted.dedup.bam",
        "sorted.dedup.bam.bai"])

# deepTools targets
DT_COR = expand(join(
    config["analysisdir"],
    config["reference"]["id"],
    "deeptools/heatmap.cor.genome.bw{binwidth}.pdf"),
     #suf = ["_10kb", "_peaks"])
    binwidth = "10000")
DT_PCA = expand(join(
    config["analysisdir"],
    config["reference"]["id"],
    "deeptools/PCA.genome.bw{binwidth}.pdf"),
     #suf = ["_10kb", "_peaks"])
    binwidth = "10000")
DT_PCA2 = expand(join(
    config["analysisdir"],
    config["reference"]["id"],
    "deeptools/plot.PCA.counts.all.genome.bw{binwidth}.pdf"),
    binwidth = "10000")
DT_FP = expand(join(
    config["analysisdir"],
    config["reference"]["id"],
    "deeptools/plot.fingerprint.bw{binwidth}.skipZeros.pdf"),
    binwidth = "500")
DT_PROF = expand(join(
    config["analysisdir"],
    config["reference"]["id"],
    "deeptools/computeMatrix/scoreMatrix.{what}.bw{binwidth}.mat.gz"),
    what = ["coverage", "log2", "subtract"],
    binwidth = config["deeptools"]["binsize"])
DT_ALL = DT_COR + DT_PCA + DT_PCA2 + DT_FP + DT_PROF
