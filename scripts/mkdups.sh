#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw504/seasonal_adaptation/log_files/array_jobs
#SBATCH --mail-type=END
#SBATCH --mail-user=rbrennan@geomar.de
#SBATCH --partition=cluster
#SBATCH --cpus-per-task=6
#SBATCH --mem=50G
#SBATCH --time=08:00:00
#SBATCH --job-name=markdups
#SBATCH --output=%x.%A.%a.out
#SBATCH --error=%x.%A.%a.err
#SBATCH --array=[1-17]

indir=$WORK/seasonal_adaptation/analysis/aligned
outdir=$WORK/seasonal_adaptation/analysis/aligned/mkdups

cd $indir
inbam=$(ls *.bam | sed -n $(echo $SLURM_ARRAY_TASK_ID)p)

echo $inbam

~/bin/sambamba-0.8.2 markdup --nthreads=6 --tmpdir=$TMPDIR $indir/$inbam $outdir/$inbam
