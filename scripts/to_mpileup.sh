#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw504/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=rbrennan@geomar.de
#SBATCH --partition=cluster
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=18:00:00
#SBATCH --job-name=to_mpileup
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err

module load samtools/1.10

cd $WORK/seasonal_adaptation/analysis/variants

samtools mpileup -q 15 -Q 0 -d 8000 -R -A -B \
	-f $WORK/seasonal_adaptation/genome/Eaffinis.Baltic.PseudoRef.Mar22.fasta \
	-b $WORK/seasonal_adaptation/analysis/variants/bamfiles.txt \
	-o all.mpileup
