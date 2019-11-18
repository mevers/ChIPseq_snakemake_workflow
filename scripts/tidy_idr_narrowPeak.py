# Snakemake Python script to convert files from an "IDR summarised peaks"
# format to a narrowPeak (BED6+4) format.
#
# Author: Maurits Evers (maurits.evers@anu.edu.au)
#
# The "IDR summarised peaks" format has the following columns
# https://github.com/nboley/idr
#   1. Chromosome
#   2. Start position (0-based)
#   3. End position (0-based, not included)
#   4. Name
#   5. Scaled IDR value min(int(log2(-125IDR), 1000); peaks with an IDR of 0
#      have a score of 1000, idr 0.05 have a score of int(-125log2(0.05)) = 540,
#      and idr 1.0 has a score of 0
#   6. Strand
#   7. Fold enrichment in merged peak region
#   8. p-value of merged peak region
#   9. q-value of merged peak region
#  10. Merged peak summit position relative to start position
#  11. -log10(local IDR value)
#  12. -log10(global IDR value)
#  13. Start position in replicate 1 (0-based)
#  14. End position in replicate 1 (0-based, not included)
#  15. Fold enrichment of peak in replicate 1
#  16. Summit position in replicate 1 relative to start position
#  17. Start position in replicate 2 (0-based)
#  18. End position in replicate 2 (0-based, not included)
#  19. Fold enrichment of peak in replicate 2
#  20. Summit position in replicate 2 relative to start position
#
# The narrowPeak format has the following columns
# http://genome.ucsc.edu/FAQ/FAQformat.html#format12
#   1. Chromosome
#   2. Start position (0-based)
#   3. End position (0-based, not included)
#   4. Name
#   5. Score for display; calculated as int(-10*log10qvalue)
#   6. Strand
#   7. Fold enrichment
#   8. -log10(p-value)
#   9. -log10(q-value)
#  10. Summit position relative to start position

import pandas as pd
import os, uuid, sys


# Snakemake parameters
infile = snakemake.input[0]
outfile = snakemake.output[0]
log = snakemake.log[0]
maxWidth = snakemake.params.maxWidth



def summarise_parameters(infile, outfile, log, maxWidth):

    # Summarise parameters
    sys.stdout = open(log, "w+")
    print("Input    : %s" % infile)
    print("Output   : %s" % outfile)
    print("maxWidth : %i" % maxWidth)


def reshape_narrowPeak(infile, outfile, maxWidth):

    # Open file
    df = pd.read_csv(infile, sep = "\t", header = None)

    # The first 10 columns are a standard narrowPeak
    df = df.iloc[:, 0:10]


    # Filter for peak width < maxLength
    df = df[(df.iloc[:, 2] - df.iloc[:, 1] <= maxWidth)]


   # Sort values
    df = df.sort_values(by = [df.columns[0], df.columns[1]])


    # Number peaks
    df.iloc[:, 3] = "IDR_peak_" + df.reset_index().index.astype(str)


    # Write to file
    df.to_csv(outfile, sep = "\t", index = False, header = False)


summarise_parameters(infile, outfile, log, maxWidth)
reshape_narrowPeak(infile, outfile, maxWidth)
