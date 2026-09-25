#!/bin/bash
#SBATCH --job-name=QC
#SBATCH -n 24
#SBATCH --partition=gen-mk-compute-1
#SBATCH --time=10:00:00
#SBATCH --mem=128G
#SBATCH --output=rrna.out
#SBATCH --error=rrna.err

module load bbmap/38.73
module load java_jdk/11.0.2

for i in `ls -1 *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//'`
do
bbduk.sh -Xmx128g in1=${i}_R1_001.fastq.gz in2=${i}_R2_001.fastq.gz out1=${i}_clean_R1.fastq.gz out2=${i}_clean_R2.fastq.gz ref=/data2/databases/rrna_silva/ribokmers.fa,adapters.fa tbo tpe ktrim=r k=31 refstats=$i.txt hdist=1 qtrim=rl trimq=15 minlength=50 maxns=1 threads=auto
done
