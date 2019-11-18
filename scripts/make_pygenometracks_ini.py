#!/usr/bin/env python

# Author: Maurits Evers

import pandas as pd
import re
import os.path, os


# Snakemake parameters
fn_bw = snakemake.input["bw"]
fn_peaks = snakemake.input["bed_peaks"]
fn_bed = snakemake.input["bed_annot"]
log = snakemake.log[0]
fn_ini = snakemake.output[0]
range = snakemake.params.range


# Summarise parameters
sys.stdout = open(log, "w+")
print("Working directory   : %s" % os.getcwd())
print("BigWig files        : %s" % ", ".join(fn_bw))
print("BED peak files      : %s" % ", ".join(fn_peaks))
print("BED annotation file : %s" % fn_bed)
print("INI output file     : %s" % fn_ini)
print("Range for y axis    : %s" % range)


# For pygenometracks it's best/easiest to work in absolute paths; otherwise
# paths need to be relative to the location of the ini file(s)
fn_bw = [os.path.join(os.getcwd(), fn) for fn in fn_bw]
fn_bed = os.path.join(os.getcwd(), fn_bed)


# Define colours
# These are the colours of the npg palette from `ggsci`
col = [
    "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
    "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF"]


# Construct Pandas DataFrame for peak files
df_peak = pd.DataFrame()
df_peak["peak"] = pd.Series(fn_peaks)
df_peak["basename"] = [os.path.basename(path) for path in df_peak["peak"]]
df_peak["cond"] = df_peak["basename"].str.replace("_IDRpeaks.narrowPeak", "")


# Infer experimental design from bigwig files and construct a Pandas DataFrame
# Assumptions:
#   1. BigWig files are labelled IP_{ID}_rep{rep_number}{anything}.bw
#   2. `cond` corresponds to the field {ID}
#   3.
df = pd.DataFrame()
df["bw"] = pd.Series(fn_bw)
df["basename"] = [os.path.basename(path) for path in df["bw"]]
df["cond"] = df["basename"].str.replace("(IP_|_rep\\d+.+$)", "")
# Infer replicate number and max replicates
df["rep"] = df["basename"].str.replace("(IP_.+_(?=rep)|_vs.+$)", "")
df["rep"] = pd.factorize(df["rep"])[0] + 1
df = df.merge(
    df.groupby("cond")["cond"].count().rename("n_rep").reset_index())
# Create group number per condition
df["group"] = pd.factorize(df["cond"])[0]
df["colour"] = df.apply(lambda x: col[x["group"]], axis = 1)
df["min_value"] = range[0]
df["max_value"] = range[1]
# Match peak files
df = df.merge(df_peak[["peak", "cond"]], on = "cond")


# Contruct axis track
entries = [
    "[x-axis]",
    "where = top",
    "title = Position",
    "fontsize = 8",
    "\n"]
track_axis = "\n".join(entries)


# Construct bw track
def format_track(row):
    entries = [
        "[bigwig %s rep%s]" % (row["cond"], row["rep"]),
        "file = %s" % row["bw"],
        "height = 3",
        "title = %s rep%s IP vs. pooled input" % (row["cond"], row["rep"]),
        "min_value = %i" % row["min_value"],
        "max_value = %i" % row["max_value"],
        "color = %s" % row["colour"],
        "\n"]
    if row["rep"] == row["n_rep"]:
        entries.extend([
            "[spacer]",
            "\n",
            "[narrowPeak %s]" % row["cond"],
            "file = %s" % row["peak"],
            "height = 1",
            "title = %s IDR peaks" % row["cond"],
            "show labels = no",
            "type = box",
            "color = %s" % row["colour"],
            "\n",
            "[spacer]",
            "\n"])
    return "\n".join(entries)

track_bw = df.apply(format_track, axis = 1).tolist()


# Construct annotation track
entries = [
    "[annotation]",
    "file = %s" % fn_bed,
    "height = 2",
	"gene rows = 5",
	"style = UCSC",
	"fontsize = 8"]

if re.search(r'U13369', fn_ini):
    entries.extend([
    "color = bed_rgb"])

track_bed = "\n".join(entries)


# Write to file
with open(fn_ini, "w") as fh:
    fh.write(track_axis)
    fh.write("".join(track_bw))
    fh.write(track_bed)

print("pygenometrack ini file successfully created.\n")
