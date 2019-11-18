# R script for ChIP-seq data snakemake analysis
# Author: Maurits Evers
# Note: DiffBind requires R >= 3.5.x

# Snakemake parameters
bam_IP <- snakemake@input["bam_IP"]
bam_IP <- c(
	"../alignment/IP_CX-5461_rep1.sorted.bam",
	"../alignment/IP_I-BET151_rep2.sorted.bam",
	"../alignment/IP_CX-5461_rep2.sorted.bam",
	"../alignment/IP_vehicle_rep1.sorted.bam",
	"../alignment/IP_doxorubicin_rep1.sorted.bam",
	"../alignment/IP_vehicle_rep2.sorted.bam",
	"../alignment/IP_CX-5461+I-BET151_rep1.sorted.bam",
	"../alignment/IP_doxorubicin_rep2.sorted.bam",
	"../alignment/IP_CX-5461+I-BET151_rep2.sorted.bam",
	"../alignment/IP_I-BET151_rep1.sorted.bam")
bam_pooled_control <- snakemake@input["bam_pooled_control"]
bam_pooled_control <- c(
	"../pooled_control_bam/control_CX-5461+I-BET151_merged.sorted.bam",
	"../pooled_control_bam/control_CX-5461_merged.sorted.bam",
	"../pooled_control_bam/control_doxorubicin_merged.sorted.bam",
	"../pooled_control_bam/control_I-BET151_merged.sorted.bam",
	"../pooled_control_bam/control_vehicle_merged.sorted.bam")
macs <- snakemake@input["macs2"]
macs <- c(
	"../macs2/callpeak/IP_CX-5461+I-BET151_rep1_vs_pooled_control_peaks.narrowPeak",
	"../macs2/callpeak/IP_doxorubicin_rep2_vs_pooled_control_peaks.narrowPeak",
	"../macs2/callpeak/IP_CX-5461+I-BET151_rep2_vs_pooled_control_peaks.narrowPeak",
	"../macs2/callpeak/IP_I-BET151_rep1_vs_pooled_control_peaks.narrowPeak",
	"../macs2/callpeak/IP_CX-5461_rep1_vs_pooled_control_peaks.narrowPeak",
	"../macs2/callpeak/IP_I-BET151_rep2_vs_pooled_control_peaks.narrowPeak",
	"../macs2/callpeak/IP_CX-5461_rep2_vs_pooled_control_peaks.narrowPeak",
	"../macs2/callpeak/IP_vehicle_rep1_vs_pooled_control_peaks.narrowPeak",
	"../macs2/callpeak/IP_doxorubicin_rep1_vs_pooled_control_peaks.narrowPeak",
	"../macs2/callpeak/IP_vehicle_rep2_vs_pooled_control_peaks.narrowPeak")


# Load libraries
library(DiffBind)
library(tidyverse)


# Infer experiment design from filenames
BAM2SampleID <- function(fn)
	str_replace(basename(fn), "(\\.sorted\\.bam)", "")
BAM2Treatment <- function(fn)
	str_replace_all(basename(fn), "(IP_|control_|_rep.+$|_merged.+$)", "")
BAM2Replicate <- function(fn)
	str_replace_all(basename(fn), "^.+rep(\\d+).+$", "\\1")
df <- tibble(bamReads = bam_IP) %>%
	mutate(
		SampleID = BAM2SampleID(bamReads),
		Treatment = BAM2Treatment(bamReads)) %>%
	left_join(
		tibble(bamControl = bam_pooled_control) %>%
			mutate(
				ControlID = BAM2SampleID(bamControl),
				Treatment = BAM2Treatment(bamControl)),
		by = "Treatment") %>%
	mutate(
		Replicate = as.integer(BAM2Replicate(bamReads)),
		PeakCaller = "narrow") %>%
	arrange(Treatment, Replicate) %>%
	as.data.frame()


# Set up differential binding analysis
res <- dba(sampleSheet = df) %>%
    dba.count(
        bUseSummarizeOverlaps = TRUE,
        bRemoveDuplicates = TRUE,
        bScaleControl = TRUE)


# Contrasts
