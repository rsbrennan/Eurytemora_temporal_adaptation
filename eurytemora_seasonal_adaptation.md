# eurytemora season adaptation

Raw data: `/gxfs_work1/geomar/smomw504/seasonal_adaptation/raw_data/`

Find sample ids that link fastq file names to sample names here: `/gxfs_work1/geomar/smomw504/seasonal_adaptation/raw_data/sample_names.txt`

lat long of samples: 54.34289 N, 9.96995 E

| Year | Sample_num | Date       | ID 		 | sample_well |
|------|------------|------------|---------------|-------------|
| 2007 | 1          | 2007-04-30 | EA_2007_T1    | B02         |
| 2007 | 2          | 2007-05-07 | EA_2007_T2    | G01         |
| 2009 | 1          | 2009-04-07 | EA_2009_T1    | A03         |
| 2009 | 2          | 2009-04-21 | EA_2009_T2    | E02         |
| 2009 | 3          | 2009-05-05 | EA_2009_T3    | A01         |
| 2009 | 4          | 2009-05-19 | EA_2009_T4    | D01         |
| 2011 | 1          | 2011-04-05 | EA_2011_T1    | F02         |
| 2011 | 2          | 2011-04-19 | EA_2011_T2    | A02         |
| 2011 | 3          | 2011-05-03 | EA_2011_T3    | G02         |
| 2015 | 1          | 2015-04-08 | EA_2015_T3    | B01         |
| 2015 | 2          | 2015-04-23 | EA_2015_T2    | H01         |
| 2015 | 3          | 2015-05-12 | EA_2015_T1    | C02         |
| 2015 | 4          | 2015-06-02 | EA_2015_T4    | F01         |
| 2022 | 1          | 2022-04-13 | EA_2022_T1    | H02         |
| 2022 | 2          | 2022-04-27 | EA_2022_T2    | E01         |
| 2022 | 3          | 2022-05-11 | EA_2022_T3    | D02         |
| 2022 | 4          | 2022-05-24 | EA_2022_T4    | C01         |

# data wrangling

data returned jan 20 2023. Sequenced at CCGA at CAU. 

library Prep was done with the Illumina DNA Prep kit: https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/illumina-dna-prep.html

Item 20060060

`We used 100ng DNA as input material and proceeeded according to the manufacturers manual.
Resulting libraries had an average fragment size of around 600 bp.
17 libraries were sequenced on 1 Lane of the NovaSeq 6000 S4 flowcell with 150bp paired-end sequencing.`


Check md5sums to make sure all downloaded ok. 

```bash

srun --pty --nodes=1 --cpus-per-task=1 --mem=1000 --time=00:40:00 /bin/bash


md5sum -c *.md5 > md5check.txt

```

All good.


## upload to ncbi:

Done on June 14, 2023

PRJNA982701

```bash

module load aspera-cli/3.9.6

ascp -i ~/aspera.openssh -QT -l100m -k1 -N '*fastq.gz$' --overwrite diff -d /gxfs_work1/geomar/smomw504/seasonal_adaptation/raw_data subasp@upload.ncbi.nlm.nih.gov:uploads/reid.brennan_gmail.com_gOcpXpBD/eurytemora2/

```

## data location:

NCBI PRJNA982701

$CEPH: `/nfs/ceph_geomar/smomw504/eurytemora_seasonal_adaptation`
$WORK: `/gxfs_work1/geomar/smomw504/seasonal_adaptation/raw_data`

## quality control and bias

fastqc on raw reads


get all programs running correctly

```bash
module load miniconda3
conda create --name trim_galore

conda activate trim_galore
conda install -c bioconda cutadapt
conda install -c bioconda fastqc
conda install -c bioconda multiqc

# Check that cutadapt is installed
cutadapt --version
# Check that FastQC is installed
fastqc -v
# Install Trim Galore
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz

```


`sbatch fastqc.sh`


look for bias in GC content with year:

```r
library(tidyverse)
dat <- read.table("~/Documents/GEOMAR/seasonal_adaptation/multiqc_general_stats.txt", header=T)
sampNm <- read.table("~/Documents/GEOMAR/seasonal_adaptation/sample_names.txt", header=F)

dat$well <- str_split_fixed(str_split_fixed(dat$Sample, "[_]", 6)[,1],"[-]", 2)[,2]

dat <- merge(x=dat, y=sampNm, by.x="well", by.y="V1")
dat$year <- substr(dat$V2, 4,7)

plot(x=dat$year, y=dat$FastQC_mqc.generalstats.fastqc.percent_gc)


## post trim
library(tidyverse)
dat <- read.table("~/Documents/GEOMAR/seasonal_adaptation/multiqc_general_stats.txt", header=T)
sampNm <- read.table("~/Documents/GEOMAR/seasonal_adaptation/sample_names.txt", header=F)

dat$well <- str_split_fixed(str_split_fixed(dat$Sample, "[_]", 6)[,1],"[-]", 2)[,2]

dat <- merge(x=dat, y=sampNm, by.x="well", by.y="V1")
dat$year <- substr(dat$V2, 4,7)

plot(x=dat$year, y=dat$FastQC_mqc.generalstats.fastqc.percent_gc)

write.table(dat, file="~/Documents/GEOMAR/seasonal_adaptation/stat_summary.txt", 
		sep="\t", quote=F, row.names=F)

```

This bias is pretty typical in older samples. We will have to deal with this- masking those sites? Pinsky's papers have dealt with these issues.


## trim

`trim.sh`


run multiqc

`multiqc.sh`


## Map

trimmed files are here: `$WORK/seasonal_adaptation/raw_data/trimmed`

There's not a super good  eurytemora genome available right now. Use the pseudoreference genome from Carol Lee's Nat Comms paper: https://www.nature.com/articles/s41467-022-31622-8 found here: https://datadryad.org/stash/dataset/doi:10.5061/dryad.r7sqv9sdz

Chromosome level assembly should be available sometime soon... at least this is what carol said in early June 2023.

I think steps are- align to pseudo ref, call snps, then use the code below to get genome position from the blast table.

https://github.com/TheDBStern/Baltic_Lab_Wild/blob/master/snp_calling/get_SNP_position_in_genome.py


Script to convert SNP positions called in one reference genome to approximate position in another genome based on blast results
parser.add_argument('-i', dest = 'snps', type = str, required=True,  help = 'input snp file, format `transcript position`')
parser.add_argument('-b', dest = 'blast', type = str, required=True,  help = 'input blast output file, outfmt 6')
parser.add_argument('-o', dest = 'out', type = str, required=True,  help = 'output converted SNP file')


### map

index genome before:

`~/bin/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index Eaffinis.Baltic.PseudoRef.Mar22.fasta.gz`

align: 

`align.sh`


### call snps

#### make mpileup

`to_mpileup.sh`

#### run varscan

`snpCall.sh`

returns 


#### filter snps


```bash
sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' all.vcf > all.4.2.vcf

sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' all.vcf > all.4.2.vcf
bcftools reheader --samples ~/seasonal_adaptation/scripts/sampleID.txt -o all.corrected.vcf all.4.2.vcf

vcftools --vcf all.corrected.vcf --minDP 30 --max-missing 1.0 --recode --recode-INFO-all --out all.nomiss

# After filtering, kept 39704 out of a possible 13150549 Sites
## why? probably bc of bad samples

conda activate tabix

zcat all.withmiss.vcftools.vcf.gz | bgzip > all.withmiss.vcftools2.vcf.gz

tabix all.withmiss.vcftools.vcf.gz 

vcftools --gzvcf all.corrected.vcf --minDP 30  --out all.withmiss --recode --recode-INFO-all

mv all.withmiss.recode.vcf all.withmiss.vcftools.vcf.gz

bcftools filter -S . -e 'FMT/DP<30 | FMT/DP>1000' -O z -o all.withmiss.vcf.gz all.corrected.vcf.gz

# from https://darencard.net/blog/2017-01-13-missing-data-proportions-vcf/
paste \
<(bcftools query -f '[%SAMPLE\t]\n' all.withmiss.vcftools.vcf.gz | head -1 | tr '\t' '\n') \
<(bcftools query -f '[%GT\t]\n' all.withmiss.vcftools.vcf.gz| awk -v OFS="\t" '{for (i=1;i<=NF;i++) if ($i == "./.") sum[i]+=1 } END {for (i in sum) print i, sum[i] / NR }' | sort -k1,1n | cut -f 2)

paste \
<(bcftools query -f '[%SAMPLE\t]\n' all.withmiss.vcf.gz | head -1 | tr '\t' '\n') \
<(bcftools query -f '[%GT\t]\n' all.withmiss.vcf.gz| awk -v OFS="\t" '{for (i=1;i<=NF;i++) if ($i == "./.") sum[i]+=1 } END {for (i in sum) print i, sum[i] / NR }' | sort -k1,1n | cut -f 2)


```
this is because the 2007 samples are terrible

% missing data < 30 or > 1000:
EA_2007_T1	0.993232
EA_2007_T2	0.993652
EA_2009_T1	0.687089
EA_2009_T2	0.45149
EA_2009_T3	0.295581
EA_2009_T4	0.158694
EA_2011_T1	0.212709
EA_2011_T2	0.154268
EA_2011_T3	0.881487
EA_2015_T1	0.142888
EA_2015_T2	0.197479
EA_2015_T3	0.533156
EA_2015_T4	0.135247
EA_2022_T1	0.140526
EA_2022_T2	0.154016
EA_2022_T3	0.13209
EA_2022_T4	0.125558

Redo filtering, but drop 2007 samples, which seemed like they totally failed.

consider EA_2011_T3 later, if still low numbers.;

```bash

vcftools --gzvcf all.withmiss.vcftools.vcf.gz --remove-indv EA_2007_T1 --remove-indv EA_2007_T2 --out all.subset --max-missing 1.0 --recode --recode-INFO-all

# After filtering, kept 1575989 out of a possible 13150549 Sites

vcftools --vcf all.subset.recode.vcf --geno-depth --out all.subset

```

get some R libraries prepped:

```bash
mkdir R_libs 
export R_LIBS=$HOME/R_libs:$R_LIBS

module load R/4.0.2 gcc/10.2.0
```

```R
library(data.table)

dat <- fread("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/all.subset.gdepth", header=T)
dat <- as.data.frame(dat)
hist(as.vector(t(dat[,3:ncol(dat)])))

quantile(as.vector(t(dat[,3:ncol(dat)])), c(0.97))

#5297


```

```bash
vcftools --vcf all.subset.recode.vcf --out all.subset --maxDP 5297  --max-missing 1.0 --recode --recode-INFO-all --min-alleles 2 --max-alleles 2 --remove-indels --out all.final

# After filtering, kept 1064850 out of a possible 1575989 Sites


conda activate tabix

bgzip all.final.recode.vcf > all.final.recode.vcf.gz

tabix all.final.recode.vcf.gz

```

#### convert to sync file for downstream stuff

```bash
module load python/3.8.6
mkdir $HOME/my_python3_env
python3 -m venv $HOME/my_python3_env/py3

module load python/3.8.6
source $HOME/my_python3_env/py3/bin/activate
module load gcc/10.2.0

pip3 install argparse

python ~/seasonal_adaptation/scripts/vcf_to_sync.py --input all.final.recode.vcf.gz --output all.final.sync

zcat all.final.recode.vcf.gz | grep TRINITY_DN89_c0_g1_i1
```
#### Now do actual filtering for af, etc. 



```R
library(poolfstat)

inNames <- read.table("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/frequency.csv", header=T, sep="\t", nrow=1)
inNames <- colnames(inNames)[grep("FREQ",colnames(inNames))]
inNames <- inNames[1:length(inNames)-1]

#inNames <- inNames[samp_drop]

dat <- vcf2pooldata(vcf.file="/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/all.final.recode.vcf",poolsizes=c(rep(100,length(inNames))),
			poolnames=inNames,min.cov.per.pool = 30,max.cov.per.pool=800,min.maf=0.01,nlines.per.readblock=1000000)

pooldata2genobaypass(dat,writing.dir="/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/")


```
filter vcf:

using script from david stern nat comms paper
```bash
python ~/seasonal_adaptation/scripts/filter_sync_by_snplist.py -i /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/all.final.sync \
	-snps /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/snpdet \
	-o /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/filtered.sync

```

### convert coordinates to full genome:

parser = argparse.ArgumentParser(description='Script to convert SNP positions called in one reference genome to approximate position in another genome based on blast results')
parser.add_argument('-i', dest = 'snps', type = str, required=True,  help = 'input snp file, format `transcript position`')
parser.add_argument('-b', dest = 'blast', type = str, required=True,  help = 'input blast output file, outfmt 6')
parser.add_argument('-o', dest = 'out', type = str, required=True,  help = 'output converted SNP file')

```bash

cut -f 1-2 /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/filtered.sync > /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/snp_coords.txt

python get_SNP_position_in_genome.py -i /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/snp_coords.txt -b /gxfs_work1/geomar/smomw504/seasonal_adaptation/genome/Eaffinis.Baltic.PseudoRef.blastn.besthits.outfmt6 -o /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/snp_positions_genome.txt


# change sync to new genome coords

cut -f 3- /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/filtered.sync | \
	paste /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/snp_positions_genome.txt - \
	> /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/filtered_genome_coords.sync

# drop the na's,  there are 5473

grep -v '^NA' /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/filtered_genome_coords.sync |\
	sort -k1,1 -k2,2n > /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/filtered_genome_coords_complete.sync

awk '!seen[$1,$2]++' /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/filtered_genome_coords_complete.sync > /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/final_genome_coords.sync


```

328,092 loci remain


now sync corresponds to affinis.Atlantic.long_read_draft.Mar22.fasta.gz


# calc Fst, allele frequencies

srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=40000 --time=03:10:00 /bin/bash


this [grenedalf](https://github.com/lczech/grenedalf) stuff maybe will no longer work. Lucas said he was changing the syntax and things in spring 2023, so probably need to make changes below.

```bash
conda activate grenedalf
module load bzip2/1.0.8

cd /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/

~/bin/grenedalf/bin/grenedalf frequency --threads 4 --log-file frequency.out \
	--sync-path /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/final_genome_coords.sync \
	--sample-name-list ~/seasonal_adaptation/scripts/sampleID_2007rm.txt \
	--write-sample-frequency \
	--write-sample-coverage \
	--write-sample-counts \
	--write-total-frequency \
	--allow-file-overwriting

### windows
~/bin/grenedalf/bin/grenedalf fst --threads 4 --log-file fst_windows.out \
	--pool-sizes 200 \
	--method unbiased-nei \
	--window-type sliding \
	--window-sliding-width 5000 \
	--window-sliding-stride 5000 \
	--omit-na-windows \
	--sync-path /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/final_genome_coords.sync  \
	--sample-name-list ~/seasonal_adaptation/scripts/sampleID_2007rm.txt \
	--allow-file-overwriting

mv fst.csv fst_windows.csv

### indiv snps
~/bin/grenedalf/bin/grenedalf fst --threads 4 --log-file fst_site.out \
	--pool-sizes 200 \
	--method unbiased-nei \
	--window-type sliding \
	--window-sliding-width 1 \
	--window-sliding-stride 0 \
	--omit-na-windows \
	--sync-path /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/final_genome_coords.sync \
	--sample-name-list ~/seasonal_adaptation/scripts/sampleID_2007rm.txt \
	--allow-file-overwriting

 mv fst.csv fst_site.csv

```


## PCA


srun --pty --x11 --nodes=1 --cpus-per-task=1 --mem=40000 --time=03:10:00 /bin/bash

module load R/4.0.2 gcc/10.2.0

fst per snp: fst_site.csv

window fst: fst_windows.csv

allele frequencies: fst_windows.csv



```R

library(data.table)

dat <- fread("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/frequency.csv", header=T, sep="\t")
dat <- as.data.frame(dat)

colin <- grep("FREQ",colnames(dat))

freqout <- dat[,colin[1:(length(colin)-1)]]

write.table(file="/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/allele_frequencies.txt", (freqout), col.names=F, row.names=F, sep="\t", quote=F)

#############################################
#############################################
##
## PCA
##
#############################################
#############################################
library(data.table)
library(ggplot2)

inNames <- fread("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/frequency.csv", header=T, sep="\t", nrow=1)
inNames <- colnames(inNames)[grep("FREQ",colnames(inNames))]
inNames <- inNames[1:length(inNames)-1]
af <- fread("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/allele_frequencies.txt", header=F, sep="\t")
af <- as.data.frame(af)
colnames(af) <- inNames

pops <- colnames(af)

pcaResult <- prcomp(t(af), scale=TRUE) 


# get proportion of total variance explained:
percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

dat.p <- data.frame(id=pops, Population=substr(pops, 4,7), 
		Treatment=substr(pops, 9,10),
		PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

a <- ggplot(dat.p, aes(PC1, PC2, fill=Treatment, shape=Population)) +
        geom_point(size=4, alpha=0.8) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
        scale_shape_manual(values=c(21,22, 23, 24))+
       ggtitle("PC1 + PC2") +
        guides(fill=guide_legend(override.aes=list(shape=21)))

ggsave("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/figures/pca.pdf",a, w=5.5, h=3.7)

##################################
# drop samples
inNames <- fread("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/frequency.csv", header=T, sep="\t", nrow=1)
inNames <- colnames(inNames)[grep("FREQ",colnames(inNames))]
inNames <- inNames[1:length(inNames)-1]

samp_drop <- grep("EA_2011_T3",inNames, invert=T )

inNames <- inNames[samp_drop]

af <- fread("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/allele_frequencies.txt", header=F, sep="\t")
af <- as.data.frame(af)
af <- af[,samp_drop]
colnames(af) <- inNames


pops <- colnames(af)

pcaResult <- prcomp(t(af), scale=TRUE) 


# get proportion of total variance explained:
percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

dat.p <- data.frame(id=pops, Population=substr(pops, 4,7), 
		Treatment=substr(pops, 9,10),
		PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

a <- ggplot(dat.p, aes(PC1, PC2, fill=Treatment, shape=Population)) +
        geom_point(size=4, alpha=0.8) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
        scale_shape_manual(values=c(21,22, 23, 24))+
       ggtitle("PC1 + PC2") +

        guides(fill=guide_legend(override.aes=list(shape=21)))

ggsave("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/figures/pca_drop1.pdf",a, w=5.5, h=3.7)


#### drop

inNames <- fread("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/frequency.csv", header=T, sep="\t", nrow=1)
inNames <- colnames(inNames)[grep("FREQ",colnames(inNames))]
inNames <- inNames[1:length(inNames)-1]

samp_drop <- grep("EA_2009_T1|EA_2011_T3",inNames, invert=T )

inNames <- inNames[samp_drop]

af <- fread("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/allele_frequencies.txt", header=F, sep="\t")
af <- as.data.frame(af)
af <- af[,samp_drop]
colnames(af) <- inNames


pops <- colnames(af)

pcaResult <- prcomp(t(af), scale=TRUE) 


# get proportion of total variance explained:
percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

dat.p <- data.frame(id=pops, Population=substr(pops, 4,7), 
		Treatment=substr(pops, 9,10),
		PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

a <- ggplot(dat.p, aes(PC1, PC2, fill=Treatment, shape=Population)) +
        geom_point(size=4, alpha=0.8) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
        scale_shape_manual(values=c(21,22, 23, 24))+
       ggtitle("PC1 + PC2") +
        guides(fill=guide_legend(override.aes=list(shape=21)))

ggsave("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/figures/pca_drop2.pdf",a, w=5.5, h=3.7)




#drop 3
inNames <- fread("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/frequency.csv", header=T, sep="\t", nrow=1)
inNames <- colnames(inNames)[grep("FREQ",colnames(inNames))]
inNames <- inNames[1:length(inNames)-1]

samp_drop <- grep("EA_2009_T1|EA_2011_T3|EA_2015_T3",inNames, invert=T )

inNames <- inNames[samp_drop]

af <- fread("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/allele_frequencies.txt", header=F, sep="\t")
af <- as.data.frame(af)
af <- af[,samp_drop]
colnames(af) <- inNames


pops <- colnames(af)
af2 <- af[ which(apply(af, 1, var) != 0), ]
nrow(af2)
af <- af2
pcaResult <- prcomp(t(af), scale=TRUE) 


# get proportion of total variance explained:
percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

dat.p <- data.frame(id=pops, Population=substr(pops, 4,7), 
		Treatment=substr(pops, 9,10),
		PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

a <- ggplot(dat.p, aes(PC1, PC2, fill=Treatment, shape=Population)) +
        geom_point(size=4, alpha=0.8) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
        scale_shape_manual(values=c(21,22, 23, 24))+
       ggtitle("PC1 + PC2") +
        guides(fill=guide_legend(override.aes=list(shape=21)))

ggsave("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/figures/pca_drop3.pdf",a, w=5.5, h=3.7)



```


There's a correlation between low quality and outliers in the plot:

EA_2009_T1
EA_2011_T3
EA_2015_T3

drop them and re-run


```R

library(data.table)
library(ggplot2)

inNames <- fread("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/frequency.csv", header=T, sep="\t", nrow=1)
inNames <- colnames(inNames)[grep("FREQ",colnames(inNames))]
inNames <- inNames[1:length(inNames)-1]

samp_drop <- grep("EA_2009_T1|EA_2011_T3|EA_2015_T3",inNames, invert=T )

inNames <- inNames[samp_drop]

af <- fread("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/allele_frequencies.txt", header=F, sep="\t")
af <- as.data.frame(af)
af <- af[,samp_drop]
colnames(af) <- inNames

pops <- colnames(af)

pcaResult <- prcomp(t(af), scale=FALSE) 
summary(pcaResult)

pc1 <- (pcaResult$rotation[,1])
pc2 <- (pcaResult$rotation[,2])
pc3 <- (pcaResult$rotation[,3])

pc_out <- data.frame(SNP= af[,1], pc1=pc1, pc2=pc2, pc3=pc3)

#write.table(file="/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/pc_loadings.txt", pc_out, quote=FALSE, row.names=FALSE, sep="\t")


# get proportion of total variance explained:
percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

dat.p <- data.frame(id=pops, Population=substr(pops, 4,7), 
		Treatment=substr(pops, 9,10),
		PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

a <- ggplot(dat.p, aes(PC1, PC2, fill=Treatment, shape=Population)) +
        geom_point(size=4, alpha=0.8) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
        scale_shape_manual(values=c(21,22, 23, 24))+
       ggtitle("PC1 + PC2") +
        guides(fill=guide_legend(override.aes=list(shape=21)))


ggsave("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/figures/pca_subset.pdf",a, w=5.5, h=3.7)


```


# below here is a mess as it currently stands.


to do next
- fst of all samples- make correlation heat map plot
- window fst, manhattan
- cvtk- not totally sure, but apply somewhow. 
- cmh, maybe just first and last
	- https://github.com/MartaPelizzola/ACER
- through time fst, 
- need to do the same type of analysis as in the drosophila seasonal adaptation stuff. 



## fst along genome

```bash
srun --pty --x11 --nodes=1 --cpus-per-task=1 --mem=20000 --time=03:10:00 /bin/bash

module load R/4.0.2 gcc/10.2.0

```

fst per snp: fst_site.csv

window fst: fst_windows.csv

```R

library(data.table)
library(tidyr)
library(corrplot)


dat <- fread("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/fst_site.csv", header=T, sep="\t")
dat <- read.csv("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/fst_windows.csv", sep="\t",header=T)

#keep only windows with reasonable number of snps

dat <- dat[(which(dat$SNPS > 5)),]

fst <- dat[,5:ncol(dat)]

fst[fst<0] <- 0
nrow(fst)
# 10832
fst <- melt(colMeans(fst, na.rm=T))
fst$id <- row.names(fst)

fst2 <- fst %>% separate_wider_delim(id, ".", names=c("id1", "id2"))

fst3 <- reshape2::acast(fst2, id1 ~ id2, value.var="value")
	fst4 <- fst3*100

pdf(file="/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/fst_windows_corr.pdf", h=8, w=8)

corrplot(fst4, method="color", type="upper", is.corr=FALSE,addCoef.col="black",
	    col=colorRampPalette(c("gray95","firebrick3"))(200),
	        tl.col="black", cl.pos ="n",tl.srt=45)
dev.off()

####

```

srun --pty --x11 --nodes=1 --cpus-per-task=1 --mem=40000 --time=03:10:00 /bin/bash

module load R/4.0.2 gcc/10.2.0


zcat affinis.Atlantic.long_read_draft.Mar22.fasta.gz > affinis.Atlantic.long_read_draft.Mar22.fasta

samtools faidx affinis.Atlantic.long_read_draft.Mar22.fasta


cut -f1-2 affinis.Atlantic.long_read_draft.Mar22.fasta.fai > affinis.Atlantic.long_read_draft.Lengths.txt

```r

library(qqman)
library(data.table)

#  SNP CHR  BP         P

dat <- read.csv("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/fst_windows.csv", sep="\t",header=T)

dat_site <- fread("/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/fst_site.csv", sep="\t",header=T)

#keep only windows with reasonable number of snps

dat <- dat[(which(dat$SNPS > 5)),]


snpCount <- table(dat$CHROM)
snpCount <- (snpCount[order(snpCount, decreasing=T)])

# drop the small ones:
snpCount <- snpCount[names(snpCount) %in% dat$CHROM]

snporder <- data.frame(CHROM = names(snpCount),
					   CHR = seq(1, length(names(snpCount)),1)
					   )


dat$BP <- dat$END
dat$SNP <- paste(dat$CHR,dat$END, sep=":")
dat$CHR <- NA

for(i in 1:nrow(snporder)){
	tmp1 <- which(dat$CHROM == snporder$CHROM[i])
	dat$CHR[tmp1] <- snporder$CHR[i]
}


library(dplyr)
dat$CHR <- as.numeric(dat$CHR)
df <- dat %>% 
   arrange(CHR,BP)


# find shared outliers:

q2009 <- quantile(df$EA_2009_T2.EA_2009_T4, 0.99)
q2011 <- quantile(df$EA_2011_T1.EA_2011_T2, 0.99)
q2015 <- quantile(df$EA_2015_T1.EA_2015_T4, 0.99)
q2022 <- quantile(df$EA_2022_T1.EA_2022_T4, 0.99)

outlier2009 <- df$SNP[which(df$EA_2009_T2.EA_2009_T4 > q2009)]
outlier2011 <- df$SNP[which(df$EA_2011_T1.EA_2011_T2 > q2011)]
outlier2015 <- df$SNP[which(df$EA_2015_T1.EA_2015_T4 > q2015)]
outlier2022 <- df$SNP[which(df$EA_2022_T1.EA_2022_T4 > q2022)]


outliers <- (intersect(c(c(outlier2009, outlier2011),outlier2015), outlier2022))

outliers <- outliers[order(outliers)]
outliers

max(df$EA_2009_T2.EA_2009_T4)
max(df$EA_2011_T1.EA_2011_T2)
max(df$EA_2015_T1.EA_2015_T4)
max(df$EA_2022_T1.EA_2022_T4)

png(file="/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/window_fst_t1_vs_t4.png", h=9, w=8, res=300, units="in")

par(mfrow=c(4,1))

manhattan(df, p = "EA_2009_T2.EA_2009_T4", logp = FALSE, ylab = "fst", genomewideline = FALSE, 
    suggestiveline = FALSE, main = "2009",
    ylim=c(0, 0.3),
    highlight=outliers)

manhattan(df, p = "EA_2011_T1.EA_2011_T2", logp = FALSE, ylab = "fst", genomewideline = FALSE, 
    suggestiveline = FALSE, main = "2011",
    ylim=c(0, 0.3),
    highlight=outliers)


manhattan(df, p = "EA_2015_T1.EA_2015_T4", logp = FALSE, ylab = "fst", genomewideline = FALSE, 
    suggestiveline = FALSE, main = "2015",
    ylim=c(0, 0.3),
    highlight=outliers)

manhattan(df, p = "EA_2022_T1.EA_2022_T4", logp = FALSE, ylab = "fst", genomewideline = FALSE, 
    suggestiveline = FALSE, main = "2022",
    ylim=c(0, 0.3),
    highlight=outliers)

dev.off()

# zoom into a few regions:

table(unlist(lapply(strsplit(as.character(outliers), ":"), '[[', 1)))
# Scz6wRH_1685;HRSCAF=1840 Scz6wRH_1693;HRSCAF=1917     Scz6wRH_23;HRSCAF=98
#                       6                       10                        6
# add CHR == #

head(df[grep("Scz6wRH_1685;HRSCAF=1840", df$CHROM),])
# chr 1
head(df[grep("Scz6wRH_1693;HRSCAF=1917", df$CHROM),])
# CHR 4
head(df[grep("Scz6wRH_23;HRSCAF=98", df$CHROM),])
# chr 3
head(df[grep("Scz6wRH_2;HRSCAF=24", df$CHROM),])
# chr 2

manhattan(subset(df, CHR == 1), p = "EA_2009_T2.EA_2009_T4", logp = FALSE, ylab = "fst", genomewideline = FALSE, 
    suggestiveline = FALSE, main = "2009",
    ylim=c(0, 0.3),
    highlight=outliers)

manhattan(subset(df, CHR == 2), p = "EA_2009_T2.EA_2009_T4", logp = FALSE, ylab = "fst", genomewideline = FALSE, 
    suggestiveline = FALSE, main = "2009",
    ylim=c(0, 0.3),
    highlight=outliers,
	xlim = c(7e7,9e7))

manhattan(subset(df, CHR == 3), p = "EA_2022_T1.EA_2022_T4", logp = FALSE, ylab = "fst", genomewideline = FALSE, 
    suggestiveline = FALSE, main = "2022",
    ylim=c(0, 0.3),
    highlight=outliers,
	xlim = c(0,2e7))

```



# covariance analysis?

want to get covariance through time for each year. 

But also look for overlap in this with the analysis across years. 


conda create --name cvtk 


conda install -c conda-forge jupyterlab 

git clone https://github.com/vsbuffalo/cvtkpy.git

python setup.py install

conda install pandas
#conda install -c anaconda numpy
conda install -c conda-forge matplotlib
conda install -c anaconda statsmodels

conda install -c conda-forge pickle5
conda install -c conda-forge regex

```bash
# on cluster
module load miniconda3
conda activate cvtk
jupyter-lab --port=6000 --no-browser

#my computer
ssh -N -L 8083:localhost:6001 smomw504@nesh-fe.rz.uni-kiel.de
http://localhost:8083/lab?token=c26c9f82253bec5f62e4ca8a063ae22a8fac55880e798082


        http://localhost:8082/lab?token=075cf34d14f5d2b73ef422734d22c150f8857c0cebfbe5b8

        http://127.0.0.1:8080/lab?token=4cc2572b48c41a90554981bf5877f73fa1dc7b97248f7210
```


make sync file with subsetted pops. equal number for all

2009: 2 and 4
2011: 1 and 2
2015: 1 and 4
2022: 1 and 4


EA_2009_T1 4
EA_2009_T2 5
EA_2009_T3 6
EA_2009_T4 7
EA_2011_T1 8
EA_2011_T2 9
EA_2011_T3 10
EA_2015_T1 11
EA_2015_T2 12
EA_2015_T3 13
EA_2015_T4 14
EA_2022_T1 15
EA_2022_T2 16 
EA_2022_T3 17
EA_2022_T4 18


```bash

cat /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/final_genome_coords.sync | cut -f 1-3,5,7,8,9,11,14,15,18 | gzip > /gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/variants_cvtk.sync.gz

/gxfs_home/geomar/smomw504/seasonal_adaptation/scripts/exp_design_cvtk_sub.csv
```

subset to just 2015 and 2022

```bash

# 2 samps frome each year
cat final_genome_coords.sync | cut -f 1-3,5,7,8,9,11,14,15,18 |gzip >variants_cvtk.sync.gz

cat final_genome_coords.sync | cut -f 1-3,5-7,12,13,14,16,17,18 | gzip > variants_2009_2015_2022.sync.gz

cat final_genome_coords.sync | cut -f 1-3,11-18 | gzip > variants_2015_2022.sync.gz

gunzip -c final_genome_coords.sync 

variants_cvtk.sync.gz

exp_design_cvtk_2009_15_22.csv


# 2 samples from each year 




```


```python

SYNC_FILE = '/gxfs_work1/geomar/smomw504/seasonal_adaptation/analysis/variants/final_genome_coords.sync'
sf = vf.SyncFile(SYNC_FILE)
sf

#
sf_gi = sf.build_gintervals()
sf_gi.seqlens = dict()
with open('/gxfs_work1/geomar/smomw504/seasonal_adaptation/genome/affinis.Atlantic.long_read_draft.Lengths.txt') as f:
    for line in f:
        seqid, length = line.strip().split('\t')
        sf_gi.seqlens[seqid] = int(length)

sf_gi

## study design

RAW_DESIGN_FILE = '/gxfs_home/geomar/smomw504/seasonal_adaptation/scripts/exp_design_cvtk_sub.csv'
sample_names = pd.read_csv(RAW_DESIGN_FILE)
#sample_re = re.compile(r"(?P<species>Dsim)_(?P<pop>Fl)_(?P<selection>Base|Hot)_(F(?P<gen>\d+)_)?(?P<rep>\d+)")
#sample_info = [{'name':d, **re.match(sample_re, d).groupdict()} for d in sample_names]
design = sample_names

design

# Our TemporalFreqs() and TiledTemporalFreqs() objects take a list of tuples (replicate, timepoint), which we create via the design DataFrame:

samples = design[['year', 'time']].copy().values.tolist()
samples


gi = sf.build_gintervals()
gi
tile_width = 1e4
tile_width_label = '1e4'
tiles = GenomicIntervals.from_tiles(gi.seqlens, width=tile_width, drop_last_tile=False)
d = TiledTemporalFreqs(tiles, freqs=sf.freqs.T, depths=sf.N.T, diploids=3000, gintervals=sf_gi, samples=samples,
                       share_first=True)



fixed = ((sf.freqs == 0.) | (sf.freqs == 1.))
# what percent of sites are fixed/lost in each replicate? similar prop?
plt.plot(fixed.mean(axis=0))

nloci = np.array([len(x) for x in d.tile_indices])

d_noempty = [ele for ele in d.tile_indices if ele != []]
print(len(d_noempty))
print(len(d.tile_indices))

print(f"mean number of loci: {nloci.mean()}")
print(f"median number of loci: {np.median(nloci)}")

```



cat ~/seasonal_adaptation/scripts/sampleID_2007rm.txt 

name,species,pop,year,time,id
1 EA_2009_T1,eurytemora,canal,2009,1,2009_1
2 EA_2009_T2,eurytemora,canal,2009,2,2009_2
3 EA_2009_T3,eurytemora,canal,2009,3,2009_3
4 EA_2009_T4,eurytemora,canal,2009,4,2009_4
5 EA_2011_T1,eurytemora,canal,2011,1,2011_1
6 EA_2011_T2,eurytemora,canal,2011,2,2011_2
7 EA_2011_T3,eurytemora,canal,2011,3,2011_3
8 EA_2015_T1,eurytemora,canal,2015,1,2015_1
9 EA_2015_T2,eurytemora,canal,2015,2,2015_2
10 EA_2015_T3,eurytemora,canal,2015,3,2015_3
11 EA_2015_T4,eurytemora,canal,2015,4,2015_4
12 EA_2022_T1,eurytemora,canal,2022,1,2022_1
13 EA_2022_T2,eurytemora,canal,2022,2,2022_2
14 EA_2022_T3,eurytemora,canal,2022,3,2022_3
15 EA_2022_T4,eurytemora,canal,2022,4,2022_4


## CMH

want to downsample so depths are biased. down to 30x for everything I think.

```bash


perl ~/bin/popoolation2/subsample-synchronized.pl --input final_genome_coords.sync --method fraction --max-coverage 1000 --target-coverage 30 --output final_genome_coords_30x.sync


# population order above. But want to use the ones that aren't funny.
# 2009 2 vs 4
# 2011 1 vs 2
# 2015 1 vs 4
# 2022 1 vs 4

perl ~/bin/popoolation2/cmh-test.pl --input final_genome_coords_30x.sync --min-coverage 30 --max-coverage 500 --min-count 1 --population 2-4,5-6,8-11,12-15 --remove-temp --output cmh_test.txt

cut -f 1-3,19- cmh_test.txt > cmh_pvals.txt

```

the last column of cmh_test.txt is the p-value.

note that some values drop out, presumably bc they don't pass filters? on 299952 lines in output.

```r
library(qqman)
#library("DataCombine")

dat <- read.table("cmh_test.txt", header=F)
frequency <- read.table("frequency.csv", sep="\t", header=T)

freq2 <- frequency[,grep("FREQ",colnames(frequency))]
freq3 <- cbind(frequency[,1:3], freq2)
colnames(dat) <- gsub(".FREQ","",colnames(freq3))

freq3$SNP <- paste(freq3$CHROM, freq3$POS, sep=":")
dat$SNP <- paste(dat$CHROM, dat$POS, sep=":")

all <- merge(freq3, dat, by="SNP", all=T)



# prep for manhattan

all$BP <- all$POS.x
all$CHR <- NA

# order the chromosomes
snpCount <- table(all$CHROM.x)
snpCount <- (snpCount[order(snpCount, decreasing=T)])


snporder <- data.frame(CHROM = names(snpCount),
					   CHR = seq(1, length(names(snpCount)),1)
					   )

for(i in 1:nrow(snporder)){
	tmp1 <- which(all$CHROM.x == snporder$CHROM[i])
	all$CHR[tmp1] <- snporder$CHR[i]
}


library(dplyr)
all$CHR <- as.numeric(dat$CHR)
df <- all %>% 
   arrange(CHR,BP)
df$TOTAL[is.na(df$TOTAL)] <- 1

# subset(df, CHR == 3)
manhattan(df, p = "TOTAL", logp = TRUE, ylab = "pvalue", genomewideline = FALSE, 
    suggestiveline = FALSE, main = "cmh")
   # ylim=c(0, 0.3),
    #highlight=outliers,


# plot loci with sig cmh with their allele frequency form t1-4 for all years.

outliers <- df[which(df$TOTAL < 0.0005),]
nrow(outliers)

data.frame(CHR = out)
# 2009 2 vs 4
# 2011 1 vs 2
# 2015 1 vs 4
# 2022 1 vs 4
t2009 <- data.frame(
			CHR = outliers$CHR,
			BP = outliers$BP,
			SNP=outliers$SNP,
			snp_id = paste(outliers$SNP, "2009", sep=":"),
			allele_frequency = c(outliers$EA_2009_T2.FREQ, outliers$EA_2009_T4.FREQ),
			time = c(rep("t1", nrow(outliers)),rep("t2", nrow(outliers))),
			year = "2009")

t2011 <- data.frame(
			CHR = outliers$CHR,
			BP = outliers$BP,
			SNP=outliers$SNP,
			snp_id = paste(outliers$SNP, "2011", sep=":"),
			allele_frequency = c(outliers$EA_2011_T1.FREQ, outliers$EA_2011_T2.FREQ),
			time = c(rep("t1", nrow(outliers)),rep("t2", nrow(outliers))),
			year = "2011")

t2015 <- data.frame(
			CHR = outliers$CHR,
			BP = outliers$BP,
			SNP=outliers$SNP,
			snp_id = paste(outliers$SNP, "2015", sep=":"),
			allele_frequency = c(outliers$EA_2015_T1.FREQ, outliers$EA_2015_T4.FREQ),
			time = c(rep("t1", nrow(outliers)),rep("t2", nrow(outliers))),
			year = "2015")

t2022 <- data.frame(
			CHR = outliers$CHR,
			BP = outliers$BP,
			SNP=outliers$SNP,
			snp_id = paste(outliers$SNP, "2022", sep=":"),
			allele_frequency = c(outliers$EA_2022_T1.FREQ, outliers$EA_2022_T4.FREQ),
			time = c(rep("t1", nrow(outliers)),rep("t2", nrow(outliers))),
			year = "2022")

allplot <- rbind(t2009, t2011,t2015,t2022)

library(ggplot2)


a <- ggplot(allplot, aes(time, allele_frequency, fill=year,color=year, shape=year, group=snp_id)) +
        geom_point(size=4, alpha=0.8, color="black") +
        geom_line() +
        #xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        #ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
        scale_shape_manual(values=c(21,22, 23, 24)) +
        facet_wrap(~SNP, scales="free_y")


```

take windowed means of cmh values.

below is from epigenetics.

need the fai file from the genome to do below. 



also calc fst means from 2009 to 2022. them compare the two.

```bash
#######################
# fixed window size:
#
#######################


# make bed file
#cat ~/Documents/GEOMAR/seasonal_adaptation/analysis/mean_fst.bed| sort -k1,1 -k2,2n >~/Documents/GEOMAR/seasonal_adaptation/analysis/mean_fst.sort.bed

cat ~/Documents/GEOMAR/seasonal_adaptation/analysis/cmh_test.txt  | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $19}'|sort -k1,1 -k2,2n  > cmh.bed


# get fst for all 4 years, plus 2009 vs 2022
# 2009 2 vs 4
# 2011 1 vs 2
# 2015 1 vs 4
# 2022 1 vs 4

cat ~/Documents/GEOMAR/seasonal_adaptation/analysis/fst_site.csv | head -n 1 | awk -v b="EA_2009_T2.EA_2009_T4" '{for (i=1;i<=NF;i++) { if ($i == b) { print i } }}'
# 20
cat ~/Documents/GEOMAR/seasonal_adaptation/analysis/fst_site.csv | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $20}' | tail -n +2 > fst_2009.bed

# 2011
cat ~/Documents/GEOMAR/seasonal_adaptation/analysis/fst_site.csv | head -n 1 | awk -v b="EA_2011_T1.EA_2011_T2" '{for (i=1;i<=NF;i++) { if ($i == b) { print i } }}'
#55
cat ~/Documents/GEOMAR/seasonal_adaptation/analysis/fst_site.csv | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $55}' | tail -n +2 > fst_2011.bed

#2015
cat ~/Documents/GEOMAR/seasonal_adaptation/analysis/fst_site.csv | head -n 1 | awk -v b="EA_2015_T1.EA_2015_T4" '{for (i=1;i<=NF;i++) { if ($i == b) { print i } }}'
# 84



cat ~/Documents/GEOMAR/seasonal_adaptation/analysis/fst_site.csv | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $20 "\t" $55 "\t" $84 "\t" $106 "\t" $54}' | tail -n +2 > fst_site_all.bed

#cat ~/tonsa_epigenetics/analysis/diff_methylation/methylation_summary.txt | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $19}' |sort -k1,1 -k2,2n > ~/Documents/GEOMAR/seasonal_adaptation/analysis/cmh.bed



# make a loop to cycle through a couple of these

win_size=2500

	bedtools makewindows -g <(cut -f 1,2 Eaffinis.Atlantic.long_read_draft.Mar22.fasta.fai) \
	-w $win_size \
	-s $win_size  \
	-i srcwinnum | sort -k1,1 -k2,2n  > ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_${win_size}kb.bed

# Join Fst values and the 'windows.bed' file
	bedtools intersect \
	-a ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_${win_size}kb.bed \
 	-b ~/Documents/GEOMAR/seasonal_adaptation/analysis/fst_site_all.bed -wa -wb | sort -k4,4  \
	> ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_fst.tab

	# cmh
	# Join Fst values and the 'windows.bed' file
	bedtools intersect \
	-a ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_${win_size}kb.bed \
 	-b ~/Documents/GEOMAR/seasonal_adaptation/analysis/cmh.bed -wa -wb | sort -k4,4  \
	> ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_cmh.tab

# Run bedtools groupby command to obtain average values of Fst for each of the windows above
#replace NA with 0

sed 's/NA/0/g' ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_fst.tab > ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_fst_narm.tab 

	bedtools groupby -i ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_fst_narm.tab \
	-g 1,2,3 \
	-c 8,9,10,11,12 \
	-o mean  > ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_mean_fst.txt

# cmh
sed 's/NA/0/g'  ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_cmh.tab | sed 's/NaN/1/g' > ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_cmh_narm.tab 

	bedtools groupby -i ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_cmh_narm.tab \
	-g 1,2,3 \
	-c 8 \
	-o mean  > ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_mean_cmh.txt

	bedtools groupby -i ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_cmh_narm.tab \
	-g 1,2,3 \
	-c 8 -o count > ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_count_cmh.txt

```

```r
#####
# get allele frequencies in these windows, just return values. 

cat ~/Documents/GEOMAR/seasonal_adaptation/analysis/fst_site.csv | head -n 1 | awk -v b="FREQ" '{for (i=1;i<=NF;i++) { if ($i == b) { print i } }}'

cat ~/Documents/GEOMAR/seasonal_adaptation/analysis/frequency.csv | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $20}' | tail -n +2 > fst_2009.bed


library(qqman)
#library("DataCombine")

frequency <- read.table("frequency.csv", sep="\t", header=T)

freq2 <- frequency[,grep("FREQ",colnames(frequency))]
freq3 <- cbind(frequency[,1:3], freq2)


bed.out <- cbind(data.frame(
		CHR= freq3$CHROM,
		START = as.character(freq3$POS-1),
		STOP = as.character(freq3$POS)),
		freq3[,(4:ncol(freq3)-1)]
		)


# write bedfile:
write.table(file="~/Documents/GEOMAR/seasonal_adaptation/analysis/frequency.bed", 
	bed.out,
	col.names=F, row.names=F, quote=F,
	sep="\t")


write.table(file="~/Documents/GEOMAR/seasonal_adaptation/analysis/frequency_header.bed", 
	colnames(bed.out),
	col.names=F, row.names=F, quote=F,
	sep="\t")

```


find intersection between af and cmh windows
```bash

# sort bed:

cat ~/Documents/GEOMAR/seasonal_adaptation/analysis/frequency.bed | sort -k1,1 -k2,2n > ~/Documents/GEOMAR/seasonal_adaptation/analysis/frequency2.bed


bedtools intersect \
	-a ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_2500kb.bed \
 	-b ~/Documents/GEOMAR/seasonal_adaptation/analysis/frequency2.bed -wa -wb | sort -k4,4  \
	> ~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_frequency.tab


```



### plot results

```r
library(qqman)
#library("DataCombine")


cmh <- read.table("~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_mean_cmh.txt", header=F)
colnames(cmh) <- c("Chrom", "start", "BP","cmh")
cmh_count <- read.table("~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_count_cmh.txt", header=F)
colnames(cmh_count) <- c("Chrom", "start", "BP","count")

nrow(cmh)
nrow(cmh_count)
cmh$count <- cmh_count$count

cmh <- cmh[cmh$count > 5,]
nrow(cmh)
# 11045

cmh$cmh[is.na(cmh$cmh)] <- 1
dat <- read.table("~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_mean_fst.txt", header=F)
colnames(dat) <- c("Chrom", "start", "BP", "fst_2009", "fst_2011", "fst_2015", "fst_2022", "fst_2009vs2022")

cmh$SNP <- paste(cmh$Chrom, cmh$BP, sep=":")
dat$SNP <- paste(dat$Chrom, dat$BP, sep=":")

all <- merge(cmh, dat, by="SNP")
# prep for manhattan

all$CHR <- NA

# order the chromosomes
snpCount <- table(all$Chrom.x)
snpCount <- (snpCount[order(snpCount, decreasing=T)])


snporder <- data.frame(CHROM = names(snpCount),
	CHR = seq(1, length(names(snpCount)),1)
					   )

for(i in 1:nrow(snporder)){
	tmp1 <- which(all$Chrom.x == snporder$CHROM[i])
	all$CHR[tmp1] <- snporder$CHR[i]
}


library(dplyr)
all$CHR <- as.numeric(all$CHR)
all$BP <- all$BP.x
df <- all %>% 
   arrange(CHR,BP)


# subset(df, CHR == 3)


quantile(df$cmh, c(0.025, 0.01))
nrow(df[which(df$cmh < 0.7311469),])

outliers <- df[which(df$cmh < 0.7311469),]
# find windows of outlier snps
consecutive_win <- as.data.frame(matrix(nrow=0, ncol=ncol(outliers)))
colnames(consecutive_win) <- colnames(outliers)

for(i in 1:(nrow(outliers)-1)){
	tmp1 <- which(df$SNP == outliers$SNP[i])
	tmp2 <- which(df$SNP == outliers$SNP[i+1])
	if(tmp2 - tmp1 ==1){
		long_win1 <- outliers[i,]
		long_win2 <- outliers[i+1,]
		consecutive_win <- rbind(consecutive_win,long_win1)
		consecutive_win <- rbind(consecutive_win,long_win2)
	}
}


consecutive_win_unique <- consecutive_win[!duplicated(consecutive_win), ]

manhattan(df, p = "cmh", logp = TRUE, ylab = "pvalue", genomewideline = FALSE, 
    suggestiveline = FALSE, main = "cmh",
    highlight=consecutive_win_unique$SNP)

manhattan(subset(df, CHR==3), p = "cmh", logp = TRUE, ylab = "pvalue", genomewideline = FALSE, 
    suggestiveline = FALSE, main = "cmh",
    highlight=consecutive_win_unique$SNP,
    xlim = c(1.12e8,1.16e8))

outliers[which(outliers$CHR == 3),]

hist(outliers$fst_2009)
hist(outliers$fst_2011)
hist(outliers$fst_2015)
hist(outliers$fst_2022)

df$fst_2022[df$fst_2022 < 0] <- 0
outliers$fst_2022[outliers$fst_2022 < 0] <- 0

plot(density(df$fst_2022))
lines(density(outliers$fst_2022),col="red")

hist(df$fst_2022, freq=F, breaks=50)
hist(outliers$fst_2022,freq=F,col=alpha("red", 0.1), breaks=50)

df$fst_2009vs2022[df$fst_2009vs2022 < 0] <- 0
outliers$fst_2009vs2022[outliers$fst_2009vs2022 < 0] <- 0
plot(density(df$fst_2009vs2022))
lines(density(outliers$fst_2009vs2022),col="red")


##############################################
##############################################
# get and plot allele freqs
##############################################
##############################################
library(dplyr)
af <- read.table("~/Documents/GEOMAR/seasonal_adaptation/analysis/windows_frequency.tab", sep="\t", header=F)

afnm <- read.table("frequency_header.bed")

colnames(af) <- c("CHR", "win_START", "win_stop", "win_snp", afnm$V1)
head(af)

af$window <- paste(
	paste(af$CHR,af$win_START, sep="_"),
	af$win_stop, sep=":")

consecutive_win_unique$window <- paste(
	paste(consecutive_win_unique$Chrom.y,consecutive_win_unique$start.y, sep="_"),
	consecutive_win_unique$BP.y, sep=":")

allwin <- merge(consecutive_win_unique,af, by="window")
nrow(allwin)
#

variants <- allwin

# polarize by rising allele
for( i in 1:nrow(variants)){
	af1 <- variants$EA_2009_T4.FREQ[i]- variants$EA_2009_T2.FREQ[i]
	af2 <- variants$EA_2011_T2.FREQ[i]- variants$EA_2011_T1.FREQ[i]
	af3 <- variants$EA_2015_T4.FREQ[i]- variants$EA_2015_T1.FREQ[i]
	af4 <- variants$EA_2022_T4.FREQ[i]- variants$EA_2022_T1.FREQ[i]

	tmpmean <- mean(c(af1, af2, af3, af4))
	# mean change < 1, flip poloarizaion
	if(tmpmean < 0){
		variants$EA_2009_T4.FREQ[i] <- 1- variants$EA_2009_T4.FREQ[i]
		variants$EA_2011_T2.FREQ[i] <- 1- variants$EA_2011_T2.FREQ[i]
		variants$EA_2015_T4.FREQ[i] <- 1- variants$EA_2015_T4.FREQ[i]
		variants$EA_2022_T4.FREQ[i] <- 1- variants$EA_2022_T4.FREQ[i]
		variants$EA_2009_T2.FREQ[i] <- 1- variants$EA_2009_T2.FREQ[i]
		variants$EA_2011_T1.FREQ[i] <- 1- variants$EA_2011_T1.FREQ[i]
		variants$EA_2015_T1.FREQ[i] <- 1- variants$EA_2015_T1.FREQ[i]
		variants$EA_2022_T1.FREQ[i] <- 1- variants$EA_2022_T1.FREQ[i]
	}

}

allwin <- variants
### make df for plotting
# 2009 2 vs 4
# 2011 1 vs 2
# 2015 1 vs 4
# 2022 1 vs 4
t2009 <- data.frame(
		window=allwin$window,
		CHR = allwin$CHR.y,
		BP = allwin$STOP,
		SNP=paste(allwin$CHR.y, allwin$STOP, sep=":"),
		snp_id = paste(
				paste(allwin$CHR.y, allwin$STOP, sep=":"), 
				"2009", sep=":"),
		allele_frequency = c(
			allwin$EA_2009_T2.FREQ, 
			allwin$EA_2009_T4.FREQ),
		time = c(rep("t1", nrow(allwin)),rep("t2", nrow(allwin))),
		year = "2009")

t2011 <- data.frame(
		window=allwin$window,
		CHR = allwin$CHR.y,
		BP = allwin$STOP,
		SNP=paste(allwin$CHR.y, allwin$STOP, sep=":"),
		snp_id = paste(
				paste(allwin$CHR.y, allwin$STOP, sep=":"), 
				"2011", sep=":"),
		allele_frequency = c(
			allwin$EA_2011_T1.FREQ, 
			allwin$EA_2011_T2.FREQ),
		time = c(rep("t1", nrow(allwin)),rep("t2", nrow(allwin))),
		year = "2011")

t2015 <- data.frame(
		window=allwin$window,
		CHR = allwin$CHR.y,
		BP = allwin$STOP,
		SNP=paste(allwin$CHR.y, allwin$STOP, sep=":"),
		snp_id = paste(
				paste(allwin$CHR.y, allwin$STOP, sep=":"), 
				"2015", sep=":"),
		allele_frequency = c(
			allwin$EA_2015_T1.FREQ, 
			allwin$EA_2015_T4.FREQ),
		time = c(rep("t1", nrow(allwin)),rep("t2", nrow(allwin))),
		year = "2015")

t2022 <- data.frame(
		window=allwin$window,
		CHR = allwin$CHR.y,
		BP = allwin$STOP,
		SNP=paste(allwin$CHR.y, allwin$STOP, sep=":"),
		snp_id = paste(
				paste(allwin$CHR.y, allwin$STOP, sep=":"), 
				"2022", sep=":"),
		allele_frequency = c(
			allwin$EA_2022_T1.FREQ, 
			allwin$EA_2022_T4.FREQ),
		time = c(rep("t1", nrow(allwin)),rep("t2", nrow(allwin))),
		year = "2022")
				

allplot <- rbind(t2009, t2011,t2015,t2022)

library(ggplot2)



for (i in unique(allplot$window)){


a <- ggplot(subset(allplot, window== i), aes(time, allele_frequency, fill=year,color=year, shape=year, group=snp_id)) +
        geom_line() +
        geom_point(size=3, alpha=0.5, color="black") +
        #xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        #ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        theme_bw() +
        scale_shape_manual(values=c(21,22, 23, 24)) +
        facet_wrap(~SNP) +
        ggtitle(i)
 	ggsave(paste0("~/Documents/GEOMAR/seasonal_adaptation/analysis/",i, ".pdf"),a, h=10, w=10)

}


```

ssh smomw504@nesh-fe.rz.uni-kiel.de
ssh -v smomw504@nesh-fe.rz.uni-kiel.de

conda create -n cvtk python=3

python -Xfrozen_modules=off -m jupyterlab_server


## PCAngsd


```bash

pcangsd

eval "$(conda shell.bash hook)" # needed for slurm
conda activate pcangsd

git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd
python setup.py build_ext --inplace
pip3 install -e .

```



cat all.subset.recode.vcf| grep -v "^#" | head


##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Raw Read Depth as reported by SAMtools">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Quality Read Depth of bases with Phred score >= 20">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">

GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR



from david stern

https://github.com/TheDBStern/Baltic_Lab_Wild/blob/master/snp_calling/bams2SNPs.commands.sh

## max depth per pop is based on overall 99% - 99.9% read count distribution (calculate_coverage_distribution_sync.py)
Rscript vcf2genobaypass.R

## Generate sync file using Popoolation

# Create the sync file, filtering bases with quality less than 20
java -ea -Xmx50g -jar <path_to_popoolation2>/mpileup2sync.jar --input lab.mpileup --output lab.sync --min-qual 20 --threads 8

### Filter the sync file to match the SNPs called with poolfstat
python filter_sync_by_snplist.py -i lab.sync -snps lab.snpdet -o lab.filtered.sync

## max depth per pop is based on overall 99% - 99.9% read count distribution (calculate_coverage_distribution_sync.py)
Rscript vcf2genobaypass.R






