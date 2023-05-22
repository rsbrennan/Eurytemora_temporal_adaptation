#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw504/seasonal_adaptation/log_files/array_jobs
#SBATCH --mail-type=END
#SBATCH --mail-user=rbrennan@geomar.de
#SBATCH --partition=cluster
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=25:00:00
#SBATCH --job-name=bwa
#SBATCH --output=%x.%A.%a.out
#SBATCH --error=%x.%A.%a.err
#SBATCH --array=[1-4]

module load samtools/1.10

indir=/gxfs_work1/geomar/smomw504/seasonal_adaptation/raw_data/trimmed
my_bwa=~/bin/bwa-mem2-2.2.1_x64-linux/bwa-mem2
bwagenind=$WORK/seasonal_adaptation/genome/Eaffinis.Baltic.PseudoRef.Mar22.fasta.gz
outdir=$WORK/seasonal_adaptation/analysis/aligned/

fq1=$(cat ~/seasonal_adaptation/scripts/redo_samples.txt | sed -n $(echo $SLURM_ARRAY_TASK_ID)p)
fq2=$(echo $fq1 | sed 's/R1_001_val_1.fq.gz/R2_001_val_2.fq.gz/g')
root=$(echo $fq1 | cut -f 8 -d "/" | cut -f 1 -d "_" | cut -f 2 -d "-" | grep -f - $WORK/seasonal_adaptation/raw_data/sample_names.txt | cut -f 2 )

rg=$(echo \@RG\\tID:$root\\tPL:Illumina\\tPU:x\\tLB:Lib1\\tSM:$root)
tempsort=$root.temp
outfile=$outdir/$root.bam

echo $SLURM_ARRAY_TASK_ID
echo $root
echo $fq1
echo $fq2
echo $rg
echo $tempsort
echo $outfile
echo $outdir

$my_bwa mem -t 4 -R $rg $bwagenind $fq1 $fq2 | \
samtools view -S -h -u - | \
samtools sort - -T $TMPDIR/$tempsort -O BAM -o $outfile
