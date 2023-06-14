#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw504/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=rbrennan@geomar.de
#SBATCH --partition=cluster
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --job-name=fastqc
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err

module load miniconda3
eval "$(conda shell.bash hook)" # needed for slurm
conda activate trim_galore

output=$WORK/seasonal_adaptation/analysis/fastqc/trimmed

cd $WORK/seasonal_adaptation/raw_data/trimmed

fastqc *fq.gz -t 16 -f fastq --noextract -o ${output}

cd $WORK/seasonal_adaptation/analysis/fastqc/trimmed

multiqc ./
