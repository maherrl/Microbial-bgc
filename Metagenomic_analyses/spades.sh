#!/bin/bash
#SBATCH --account=crobe
#SBATCH --partition=longfat
#SBATCH --job-name=spades
#SBATCH --output=spades_BIGC_indiv2.out
#SBATCH --error=spades_BIGC_indiv2.err
#SBATCH --time=10000
#SBATCH --mem=250000M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6

base="/home/rmaher2/crobe/NEON/2019"
in="${base}/clean"
out="${base}/assembly"

module load spades/3.13.0

mkdir ${out}/BIGC.20190709.EPILITHON_indiv2

spades.py -o ${out}/BIGC.20190709.EPILITHON_indiv2 \
    -1 ${in}/BIGC.20190709.EPILITHON.2.DNA-DNA1_clean_R1.fastq -2 ${in}/BIGC.20190709.EPILITHON.2.DNA-DNA1_clean_R2.fastq \
    -t 6 \
    --tmp-dir /tmp \
    --only-assembler
