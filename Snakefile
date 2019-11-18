# Snakemake workflow for the analysis of ChIP-seq data
#
# Author: Maurits Evers
# License: GPLv3
# Original date: 14-10-2016
# Last changed: 15-08-2017

from os.path import join
import re
import glob

def message(mes):
  sys.stderr.write("|--- " + mes + "\n")

#######################################################
################## Global variables ###################
#######################################################

# Config files
configfile: "config.yaml"

# Working directory
workdir: config["basedir"]

#######################################################
######################### Targets #####################
#######################################################

# Reference sequence
REF = join(
    config["refdir"],
    config["reference"]["id"],
    config["reference"]["filename"])

# bowtie2 index
IDX = expand(re.sub("fa", "{idx}.bt2", REF), idx = range(1,5))

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


DT_RATIO = expand(join(
    config["analysisdir"],
    config["reference"]["id"],
    "deeptools/bamCompare/{operation}/{name}_vs_pooled_control.normSES.bw10.bw"),
    operation = ["log2", "ratio", "subtract"],
    name = [sample
        for treatment in config["samples"].keys()
            for sample in config["samples"][treatment]["IP"]])

IDR = [join(
    config["analysisdir"],
    config["reference"]["id"],
    "idr",
    treatment + "_IDRsummarised.narrowPeak")
    for treatment in config["samples"].keys()]

PGT = expand(join(
    config["analysisdir"],
    config["reference"]["id"],
    "pygenometracks/{operation}_IP_vs_pooled_control_{region}_bw10.{suffix}"),
    operation = ["log2"],
    region = [
        "6:135165569-135234910",   # MYB
        "6:73509724-73527058",     # EEF1A1
        "9:97341742-97566477",     # TDRD7
        "8:127724462-127753922",   # MYC
        "U13369.1:1-42998"],       # rDNA
    suffix = ["pdf", "png"])

INI = expand(join(
    config["analysisdir"],
    config["reference"]["id"],
    "pygenometracks",
    "tracks_{region}_{operation}_bw{binwidth}.ini"),
    operation = ["ratio"],
    region = ["U13369.1", "genome"],
    binwidth = 10)


for smp in INI:
    message("Sample " + smp + " will be created")

#######################################################
###################### Includes #######################
#######################################################


include: "rules/reference.rules"
include: "rules/process_fastq.smk"
include: "rules/alignment.smk"
include: "rules/bam.smk"
include: "rules/deeptools_data_prep.smk"
include: "rules/R.rules"
include: "rules/macs2.smk"
include: "rules/idr.smk"
include: "rules/pygenometracks.smk"
#include: "rules/deeptools_plots.smk"


#######################################################
######################## Rules ########################
#######################################################

# Input fastq files
rule all:
    input:
        DT_RATIO + IDR + PGT
