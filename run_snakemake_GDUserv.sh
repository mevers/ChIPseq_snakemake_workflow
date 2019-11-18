#!/bin/bash

echo "Plotting dag's..."
snakemake --dag | dot -Tpdf > dag.pdf
snakemake --rulegraph | dot -Tpdf > rulegraph.pdf
echo "[DONE]"

echo "Running main snakemake workflow ..."
snakemake --configfile config.yaml \
	  --snakefile Snakefile \
	  --jobs 10 \
	  --printshellcmds \
      --reason \
      --cluster-config cluster.yaml \
	  --cluster "qsub -pe threads {cluster.threads} \
                          -q {cluster.queue} \
                          -l virtual_free={cluster.virtual_free} \
                          -l h_vmem={cluster.h_vmem} \
                          -o {cluster.outstream} \
                          -e {cluster.errorstream}"
echo "[DONE]"

#echo "Running snakemake generate_bowtie2_summary_stats..."
#snakemake --configfile config.yaml \
#      --snakefile Snakefile \
#      --jobs 1 \
#	  --printshellcmds \
#      --reason \
#      --cluster "qsub -pe threads {cluster.threads} \
#                          -q {cluster.queue} \
#                          -l virtual_free={cluster.virtual_free} \
#                          -l h_vmem={cluster.h_vmem} \
#                          -o {cluster.outstream} \
#                          -e {cluster.errorstream}" \
#      --cluster-config cluster.yaml \
#      generate_bowtie2_summary_stats
#echo "[DONE]"
#
#echo "Running snakemake generate_MarkDuplicates_summary_stats..."
#snakemake --configfile config.yaml \
#      --snakefile Snakefile \
#      --jobs 1 \
#	  --printshellcmds \
#      --timestamp \
#      --reason \
#      --cluster "qsub -pe threads {cluster.threads} \
#                          -q {cluster.queue} \
#                          -l virtual_free={cluster.virtual_free} \
#                          -l h_vmem={cluster.h_vmem} \
#                          -o {cluster.outstream} \
#                          -e {cluster.errorstream}" \
#      --cluster-config cluster.yaml \
#      generate_MarkDuplicates_summary_stats
#echo "[DONE]"
