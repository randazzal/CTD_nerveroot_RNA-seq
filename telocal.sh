#!/bin/bash
#SBATCH --job-name=local_NR
#SBATCH --cpus-per-task=12
#SBATCH --partition=gen-mk-compute-1
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --array=1-22%5
#SBATCH --output=logs/local_%A_%a.out
#SBATCH --error=logs/local_%A_%a.err
#conda activate telocal

name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)


TElocal -b ${name}_Aligned.out_sorted.bam \
	--GTF /data2/lackey_lab/GenomeReferences/hisat_index/gencode.v39.annotation.gtf \
	--TE /data2/lackey_lab/DownloadedSequenceData/randazza/TElocal/GRCh38_GENCODE_rmsk_TE.gtf.locInd \
	--sortByPos \
	--project TElocal/${name}_telocale
