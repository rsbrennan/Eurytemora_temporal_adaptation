#!/bin/bash
#SBATCH -D /gxfs_home/geomar/smomw504/seasonal_adaptation/log_files/
#SBATCH --mail-type=END
#SBATCH --mail-user=rbrennan@geomar.de
#SBATCH --partition=cluster
#SBATCH --cpus-per-task=1
#SBATCH --mem=28G
#SBATCH --time=18:00:00
#SBATCH --job-name=align_stats
#SBATCH --output=%x.%A.out
#SBATCH --error=%x.%A.err

module load samtools

cd $WORK/seasonal_adaptation/analysis/aligned/mkdups

echo "sample,total_reads,total_mapped,map_percent,mapped_q20,map_percent_q20,duplicates,duplicate_percent"> $WORK/seasonal_adaptation/analysis/aligned/mkdups/count.aligned.txt

for i in $(ls *.bam | cut -f -1 -d "." | uniq )

do {

  samtools flagstat ${i}.bam> tmp.txt

  TOTAL=$(grep "in total" tmp.txt | cut -f 1 -d " ")
  DUPS=$(grep "duplicates" tmp.txt | cut -f 1 -d " ")
  MAPPED=$(grep "mapped (" tmp.txt | cut -f 1 -d " ")
  PERCENT=$(grep "mapped (" tmp.txt| cut -f 2 -d "(" | cut -f 1 -d "%")
  MAPPED_q=$(samtools view -F 4 -q 20 ${i}.bam | wc -l )
  PERCENT_q=$(echo "scale=2 ; $MAPPED_q / $TOTAL" | bc)
  PERCENT_dup=$(echo "scale=2 ; $DUPS / $TOTAL" | bc)


  echo ${i},$TOTAL,$MAPPED,$PERCENT,$MAPPED_q,$PERCENT_q,$DUPS,$PERCENT_dup

} >> $WORK/seasonal_adaptation/analysis/aligned/mkdups/count.aligned.txt

done

  rm tmp.txt
