#!/bin/bash
#SBATCH --account=crobe
#SBATCH --partition=longfat
#SBATCH --job-name=bwt
#SBATCH --output=bwt3.out
#SBATCH --error=bwt3.err
#SBATCH --time=1500
#SBATCH --mem=50000M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20

NUM_THREADS=20
base="/home/rmaher2/crobe/neon/2019"
in="${base}/clean"
contigs="${base}/contigs_reformatted"
out="${base}/mapping_reps"

module load samtools/1.5
module load bowtie2/2.4.4

for sample in `awk '{print $1}' ${in}/sample_reps.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi
    
    # you need to make sure you "ls 01_QC/*QUALITY_PASSED_R1*" returns R1 files for all your samples in samples.txt
    R1s=`ls ${in}/$sample.DNA-DNA1_clean_R1.fastq | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
    R2s=`ls ${in}/$sample.DNA-DNA1_clean_R2.fastq | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))'`
    assembly=${sample::-2}

    bowtie2-build ${contigs}/${assembly}_benthic.fasta ${contigs}/${assembly}_benthic.db --threads $NUM_THREADS
    bowtie2 --threads $NUM_THREADS -x ${contigs}/$assembly/contigs.db -1 $R1s -2 $R2s --no-unal -S ${out}/$sample.sam
    samtools view -F 4 -bS ${out}/$sample.sam > ${out}/$sample-RAW.bam
    anvi-init-bam ${out}/$sample-RAW.bam -o ${out}/$sample.bam
done
