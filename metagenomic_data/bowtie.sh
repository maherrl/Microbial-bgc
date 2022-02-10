#!/bin/bash
#SBATCH --account=crobe
#SBATCH --partition=fat
#SBATCH --job-name=bwt
#SBATCH --output=bwt.out
#SBATCH --error=bwt.err
#SBATCH --time=1000
#SBATCH --mem=25000M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15

base="/home/rmaher2/crobe/neon/2019"
in="${base}/anvio/contigs"
reads="${base}/clean"
out="${base}/anvio/mapping/BLDE"

module load bowtie2/2.4.4
module load samtools/1.14

#bowtie2-build ${in}/BLDE_spades_co.fasta ${in}/bwt_dbs/BLDE_spades_co.db --threads 15

for sample in `awk '{print $1}' ${reads}/BLDE.samples.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi
    bowtie2 --threads 15 \
    -x ${in}/bwt_dbs/BLDE_spades_co.db \
    -1 ${reads}/${sample}_R1.fastq \
    -2 ${reads}/${sample}_R2.fastq \
    --no-unal \
    -S ${out}/${sample}.sam
    
    samtools view -F 4 -bS ${out}/${sample}.sam > ${out}/${sample}-RAW.bam
    anvi-init-bam ${out}/${sample}-RAW.bam -o ${out}/${sample}.bam

done
