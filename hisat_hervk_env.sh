#!/bin/bash
#SBATCH --job-name=hervk_env_hisat
#SBATCH -n 24
#SBATCH --partition=gen-mk-compute-1
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --output=hisat.out
#SBATCH --error=hisat.err

# Load software
module load hisat2/2.1.0

for i in *_R1.fastq.gz; do name=$(basename ${i} _R1.fastq.gz);
hisat2 --no-spliced-alignment -k 10 -x hervk3_index \
-1 ${i} \
-2 ${name}_R2.fastq.gz \
-S ${name}_hervk3.sam
done
