#!/bin/bash
#SBATCH -D /home/smomw504/seasonal_adaptation/log_files/array_jobs
#SBATCH --mail-type=END
#SBATCH --mail-user=rbrennan@geomar.de
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --job-name=trim
#SBATCH --output=%x.%A.%a.out
#SBATCH --error=%x.%A.%a.err
#SBATCH --array=[1-17]%6

module load miniconda3
eval "$(conda shell.bash hook)" # needed for slurm
conda activate trim_galore

# Search all the fastq files from the "indir" directory and generate the array
cd /work_beegfs/smomw504/seasonal_adaptation/jf_brennan_geomar_328_hmvl2dsx5

# identify the fq files.
fq1=$(ls *_R1_001.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
fq2=$(echo $fq1 | sed 's/_R1/_R2/g')


echo $SLURM_ARRAY_TASK_ID
echo $fq1
echo $fq2

~/bin/TrimGalore-0.6.6/trim_galore --length 40 --nextera --gzip --output_dir trimmed --paired ${fq1} ${fq2}
