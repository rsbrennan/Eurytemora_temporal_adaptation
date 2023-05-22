#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw504/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=rbrennan@geomar.de
#SBATCH --partition=cluster
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=38:00:00
#SBATCH --job-name=snpCall
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err

module load openjdk

cd $WORK/seasonal_adaptation/analysis/variants

java -jar ~/bin/varscan-2.4.5/VarScan.v2.4.5.jar mpileup2cns all.mpileup --min-coverage 20 --min-avg-qual 20 --min-reads2 4 --p-value 0.1 --min-var-freq 0.001 --output-vcf 1 --variants > all.vcf
