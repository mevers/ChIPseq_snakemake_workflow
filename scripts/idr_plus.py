#!/usr/bin/env python3

from pathlib import Path
from shutil import copy2
import os, uuid, sys
import subprocess


# Snakemake parameters
sn_in = snakemake.input
sn_out = snakemake.output[0]
sn_log = snakemake.log[0]


# Create folder
dir = os.path.join(os.path.dirname(sn_out), "idr_plus")
if not os.path.exists(dir):
    os.makedirs(dir)


def summarise_parameters(infiles, outfile, log):

    # Summarise parameters
    sys.stdout = open(log, "w+")
    print("Input    : %s" % ", ".join(infiles))
    print("Output   : %s" % outfile)


def parse_args():

    import argparse

    parser = argparse.ArgumentParser(
        formatter_class = argparse.RawTextHelpFormatter,
        description = """
        IDR (Irreproducible Discovery Rate) analysis with > 2 replicates

        This is wrapper method around `idr` that iteratively performs pairwise
        IDR analyses and collapses results into one IDR summarised peak file.

        Example 1:
            idr_plus --samples peaks1.narrowPeak peaks2.narrowPeak peaks3.narrowPeak

        This will perform the following IDR analyses:
            1. IDR(peaks1 & peaks2) giving peaks12
            2. IDR(peaks1 & peaks3) giving peaks13
            3. IDR(peaks12 & peaks13) giving peaks123

        Example 2:
            idr_plus --samples peaks1.narrowPeak peaks2.narrowPeak peaks3.narrowPeak
                peaks4.narrowPeak

        This will perform the following IDR analyses:
            1. IDR(peaks1 & peaks2) giving peaks12
            2. IDR(peaks1 & peaks3) giving peaks13
            3. IDR(peaks1 & peaks4) giving peaks14
            4. IDR(peaks12 & peaks13) giving peaks123
            5. IDR(peaks12 & peaks14) giving peaks124
            6. IDR(peaks123 & peaks124) giving peaks1234
        """
    )

    parser.add_argument(
        "--samples", "-s",
        nargs = "*",
        required = True,
        help = "Files containing peaks and score [in narrowPeak format]")

    parser.add_argument(
        "--output-file", "-o",
        nargs = 1,
        required = True,
        help = "File containing the final IDR-summarised peaks")


    args = parser.parse_args()

    return args


# Create copies of input files in temporary folder
def init_files(smpl):
    # Copy files
    fn = list()
    for i in range(len(smpl)):
        dest = os.path.join(dir, "%s.narrowPeak" % (i + 1))
        res = copy2(smpl[i], dest)
        fn.append(dest)
    print("Created temporary copies in folder %s." % dir)
    # Return list of new filenames
    return fn


# Split a string into characters
def split(s):
    return [c for c in s]


# Recursively do pairwise IDRs
def recursive_pairwise_op(x):

    lst = list()
    while len(x) > 1:

        # Generate output file
        id = sorted(set(split(Path(x[0]).stem) + split(Path(x[1]).stem)))
        out = os.path.join(dir, "%s.narrowPeak" % "".join(id))

        lst.append(out)
        cmd = """
            idr \
            --samples %s \
            --input-file-type narrowPeak \
            --output-file %s \
            --plot \
            --peak-merge-method avg \
            --verbose
            """ % (" ".join(x[0:2]), lst[-1])
        subprocess.run(cmd, shell = True)
        #subprocess.run(["echo", cmd])
        del(x[1])


    if len(lst) > 1:
        return recursive_pairwise_op(lst)
    else:
        return lst[0]


def main():

    #args = parse_args()
    #smpl = args.samples
    #out = args.output_file[0]

    summarise_parameters(sn_in, sn_out, sn_log)

    fn = init_files(sn_in)

    res = recursive_pairwise_op(fn)

    status = copy2(res, sn_out)
    print("Final IDR-summarised peak file in %s." % sn_out)


main()
