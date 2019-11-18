#!/usr/bin/env python

# Author: Maurits Evers

import sys
import os
import gzip
import re
import json
import tqdm
import struct
import itertools


# Snakemake parameters
fn1 = snakemake.input[0]
fn2 = snakemake.input[1]
dir = snakemake.params.dir


# process returns a dictionary
def process(lines = None):
    ks = ["ID", "sequence", "optional", "quality"]
    return {key: val for key, val in zip(ks, lines)}


# Test if file exists
def file_exist(fn):
    if not os.path.exists(fn):
        raise SystemError("Error: File %s does not exist\n" % fn)


# Remove substring
def remove_string(string, start, end):
    return string[:start] + string[end:]


# Check that input files exist
file_exist(fn1)
file_exist(fn2)


# Generate output files
fn_out1 = re.sub(".fastq.gz", ".polyTA_reads_only.fastq.gz", fn1)
fn_out2 = re.sub(".fastq.gz", ".polyTA_reads_only.fastq.gz", fn2)
fn_out1 = os.path.join(dir, os.path.basename(fn_out1))
fn_out2 = os.path.join(dir, os.path.basename(fn_out2))


# Read in records from both sequence files in parallel
n = 4
dict = {}
totalsize = os.path.getsize(fn1)
with tqdm.tqdm(total = totalsize) as pbar:
    with open(fn1, "rb") as fh1, open(fn2, "rb") as fh2, \
         open(fn_out1, "wb") as fh_out1, open(fn_out2, "wb") as fh_out2:

        # Input files
        in1 = gzip.GzipFile(fileobj = fh1)
        in2 = gzip.GzipFile(fileobj = fh2)

        # Output files
        out1 = gzip.GzipFile(fileobj = fh_out1)
        out2 = gzip.GzipFile(fileobj = fh_out2)


        nrec_pre = 0
        nrec_post = 0
        lines1 = []
        lines2 = []
        size0 = 0
        size1 = 0

        for line1, line2 in zip(in1, in2):
            lines1.append(line1.rstrip().decode("ascii"))
            lines2.append(line2.rstrip().decode("ascii"))
            if len(lines1) == n:
                nrec_pre += 1
                record1 = process(lines1)
                record2 = process(lines2)
                regex = re.compile("((?<=[ACG])T{3,}$|^A{3,}(?=[CGT]))")
                m1 = regex.search(record1["sequence"])
                m2 = regex.search(record2["sequence"])

                if m1 or m2:

                    if m1 and not m2:
                        nrec_post += 1

                    if m2 and not m1:
                        nrec_post += 1

                    if m1 and m2:
                        nrec_post += 1

                    out1.write((record1["ID"] + "\n").encode())
                    out1.write((record1["sequence"] + "\n").encode())
                    out1.write((record1["optional"] + "\n").encode())
                    out1.write((record1["quality"] + "\n").encode())
                    out2.write((record2["ID"] + "\n").encode())
                    out2.write((record2["sequence"] + "\n").encode())
                    out2.write((record2["optional"] + "\n").encode())
                    out2.write((record2["quality"] + "\n").encode())


                lines1 = []
                lines2 = []

            # Update status bar
            size1 = fh1.tell()
            pbar.update(size1 - size0)
            size0 = size1

        # Close files
        out1.close()
        out2.close()

print("Sequencing files : %s, %s" % (fn1, fn2))
print("Pre-trimming     : %i records per file." % nrec_pre)
print("Post-trimming    : %i records per file." % nrec_post)
print("Outout files     : %s, %s" % (fn_out1, fn_out2))
