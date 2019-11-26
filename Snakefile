# Snakemake workflow for the analysis of ChIP-seq data
#
# Author: Maurits Evers
# License: GPLv3
# Original date: 14-10-2016
# Last changed: 15-08-2017

from os.path import join
import re
import glob
import itertools

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

# DeepTools BigWig targets
DT_RATIO = expand(join(
    config["analysisdir"],
    config["reference"]["id"],
    "deeptools/bamCompare/{operation}",
    "{name}_vs_pooled_control.norm{scale_method}.bw10.bw"),
    operation = ["log2", "ratio", "subtract"],
    name = [sample
        for treatment in config["samples"].keys()
            for sample in config["samples"][treatment]["IP"]],
    scale_method = config["deeptools"]["scale_method"])

# IDR targets
IDR = [join(
    config["analysisdir"],
    config["reference"]["id"],
    "idr",
    treatment + "_" + comparison + ".narrowPeak")
    for treatment in config["samples"].keys()
        for comparison in ["_vs_".join(
            re.sub(r"IP_(\w+_)*", "IDR", x) for x in w)
            for w in itertools.combinations(
                config["samples"][treatment]["IP"], 2)]]

# PyGenomeTracks targets
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
        "BK000964.3.1:1-45306"],   # rDNA
    suffix = ["pdf", "png"])


INI = expand(join(
    config["analysisdir"],
    config["reference"]["id"],
    "pygenometracks",
    "tracks_{region}_{operation}_bw{binwidth}.ini"),
    operation = ["ratio"],
    region = ["BK000964.3", "genome"],
    binwidth = 10)


#for smp in IDR:
#    message("Sample " + smp + " will be created")

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
        DT_RATIO + IDR + INI
