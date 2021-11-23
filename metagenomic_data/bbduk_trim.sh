#!/bin/bash
#SBATCH --account=crobe
#SBATCH --partition=short
#SBATCH --job-name=bbduk
#SBATCH --output=bbduk.out
#SBATCH --error=bbduk.err
#SBATCH --time=1440
#SBATCH --mem=20000M
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5

module load easybuild
module load icc/2017.1.132-GCC-6.3.0-2.27
module load impi/2017.1.132
module spider BBMap/37.25-Java-1.8.0_162

base="/home/rmaher2/crobe/NEON/2019"
out="${base}/trim_bb"

cd ${base}

for f in *_R1.fastq

do
	n=${f%%_R1.fastq}
	bbduk.sh -Xmx1g threads=5 in1=${base}/${n}_R1.fastq in2=${base}/${n}_R2.fastq \
	out1=${out}/${n}_clean_R1.fastq out2=${out}/${n}_clean_R2.fastq \
	ref=${out}/adapters.fa ktrim=r k=12 mink=11 hdist=1 tpe tbo
done
