#!/usr/bin/env python

# Author: Maurits Evers

import pandas as pd
import os, uuid, sys


# Snakemake parameters
infile = snakemake.input[0]
outfile = snakemake.output[0]
log = snakemake.log[0]
thresholdEF = snakemake.params.thresholdEF


# Open log file
sys.stdout = open(log, "w+")
print("Input                       : %s" % infile)
print("Output                      : %s" % outfile)
print("EF threshold                : %3.2f" % thresholdEF)


# Open file
df = pd.read_csv(infile, sep = "\t", header = None)
print("Number of entries in input  : %i" % len(df.index))


# Only keep entries where the EF is > thresholdEF
df = df[(df.iloc[:, 6] >= thresholdEF)]


# Write to file
df.to_csv(outfile, sep = "\t", index = False, header = False)
print("Number of entries in output : %i" % len(df.index))
