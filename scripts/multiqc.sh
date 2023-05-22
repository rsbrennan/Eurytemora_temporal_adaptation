#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw504/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=rbrennan@geomar.de
#SBATCH --partition=cluster
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --job-name=multiqc
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err

module load miniconda3
eval "$(conda shell.bash hook)" # needed for slurm
conda activate trim_galore

cd /gxfs_work1/geomar/smomw504/seasonal_adaptation/raw_data/trimmed

multiqc ./
