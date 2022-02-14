#!/bin/bash
#SBATCH --account=crobe
#SBATCH --partition=longfat
#SBATCH --job-name=bwt
#SBATCH --output=bwt.out
#SBATCH --error=bwt.err
#SBATCH --time=10000
#SBATCH --mem=500000M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15

NUM_THREADS=15
base="/home/rmaher2/crobe/neon/2019"
in="${base}/clean"
contigs="${base}/assembly"
out="${base}/mapping"

module load bowtie2/2.4.4

for sample in `awk '{print $1}' ${in}/samples.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi
    
    # you need to make sure you "ls 01_QC/*QUALITY_PASSED_R1*" returns R1 files for all your samples in samples.txt
    R1s=`ls ${in}/$sample*.DNA-DNA1_clean_R1.fastq | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
    R2s=`ls ${in}/$sample*.DNA-DNA1_clean_R2.fastq | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
    
    bowtie2-build ${contigs}/$sample/contigs.fasta ${contigs}/$sample/contigs.db --threads $NUM_THREADS
    bowtie2 --threads $NUM_THREADS -x ${contigs}/$sample/contigs.db -1 $R1s -2 $R2s --no-unal -S ${out}/$sample.sam
    samtools view -F 4 -bS ${out}/$sample.sam > ${out}/$sample-RAW.bam
    anvi-init-bam ${out}/$sample-RAW.bam -o ${out}/$sample.bam
done
