![Image 1](bedbugs.png)

Contact : chloe.haberkorn@univ-lyon1.fr

### Table of Contents

- **[Pool-seq data processing](#Pool-seq-data-processing)**
	- [Installing tools](#Installing-tools)
	- [Getting the data](#Getting-the-data)
	- [Trimming](#Trimming)
	- [Removing duplicates](#Removing-duplicates)
	- [Mapping](#Mapping)
	- [Analysing coverage](#Analysing-coverage)

- **[Detecting Single Nucleotide Polymorphism](#Detecting-Single-Nucleotide-Polymorphism)**
	- [Installing tools](#Installing-tools)
	- [Overall SNPs analyzes](#Overall-SNPs-analyzes)

 - **[Selecting candidate SNPs](#Selecting-candidate-SNPs)**
	- [Installing tools](#Installing-tools)
	- [Overall SNPs analyzes](#Overall-SNPs-analyzes)



- **[Genetic polymorphism within populations](#Genetic-polymorphism-within-populations)**
	- [Installing tools](#Installing-tools)
	- [Compute pi and Tajima's D](#Compute-pi-and-Tajimas-D)

- **[Evidence of selection](#Evidence-of-selection)**
	- [Installing tools](#Installing-tools)
	- [Identifying genetic markers under selection](#Identifying-genetic-markers-under-selection)
	- [Mapping areas under selection](#Mapping-areas-under-selection)

## Pool-seq data processing

The goal is first to map *Cimex lectularius* PoolSeq samples (London Lab, London Field, German Lab and Sweden Field - pools of 30 individuals) on reference genome.
We used the recent reference genome and annotation, avalaible here: https://www.ncbi.nlm.nih.gov/assembly/GCF_000648675.2

We will have to download a few softs.

### Install tools

Here are the tools and versions used: 
- FastQC 
- Trimmomatic v0.39
- FastUniq v1.1
- BWA v0.74
- Samtools v1.9
- Bedtools v2.29.1

They will be store in /your-path/Tools.

Example for BWA:
``` 
git clone https://github.com/lh3/bwa.git
cd bwa
make
/your-path/Tools/bwa/bwa index # Check that the soft is working 
```

For Trimmomatic, we will also have to download the associated adapters. To do so : in fastq file, look at overrepresented sequences -> "TruSeq Adapter". These primers come from the TruSeq-3 library.

### Getting the data

Raw sequences (fastq.gz files):
```
mkdir /your-path/PoolSeq_Clec/Raw_Clec
```

Reference genome:
```
mkdir /your-path/PoolSeq_Clec/Ref_Clec
```

### Trimming

Keep default parameters, except for:
SLIDINGWINDOW:4:15 -> 20 was too high, use 15 instead

ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads 
  - TruSeq3-PE.fa = adaptaters to remove
  - 2 = seed mismatches
  - 30 = palindrome clip threshold
  - 10 = simple clip threshold
  - 2 = min Adapter Length

```
mkdir /your-path/PoolSeq_Clec/Trimmed
cd /your-path/PoolSeq_Clec/Trimmed

DIR=/your-path/PoolSeq_Clec
DIRSOFT=/your-path/Tools
DIRFASTQ="$DIR"/Raw_Clec

/usr/local/jre1.8.0_202/bin/java -jar "DIRSOFT"/trimmomatic-0.39.jar PE -phred33 -trimlog LL_trim.log \
"DIRFASTQ"/LL_R1.fastq.gz "DIRFASTQ"/LL_R2.fastq.gz \
LL_R1_paired.fq.gz LL_R1_unpaired.fq.gz LL_R2_paired.fq.gz LL_R2_unpaired.fq.gz \
ILLUMINACLIP:"DIRSOFT"/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 

gunzip LF_*_paired.fq.gz
```

Check filtering quality using fastq:
```
/your-path/Tools/FastQC/fastqc --java /usr/local/jre1.8.0_202/bin/java LF_R1_paired.fq
/your-path/Tools/FastQC/fastqc --java /usr/local/jre1.8.0_202/bin/java LF_R2_paired.fq
```

### Removing duplicates

Create a text file with both input - Example for input_LL.txt:
```
LL_R1_paired.fq
LL_R2_paired.fq
```

Using FastUniq:
```
cd /your-path/PoolSeq_Clec/Trimmed
gunzip *_paired.fq.gz

/your-path/Tools/FastUniq/source/fastuniq -i input_LL.txt -t q -o LL_dup_1.fastq -p LL_dup_2.fastq -c 1
```

Check how many duplicate sequences have been deleted:
```
/your-path/Tools/FastQC/fastqc --java /usr/local/jre1.8.0_202/bin/java LL_dup_1.fastq
/your-path/Tools/FastQC/fastqc --java /usr/local/jre1.8.0_202/bin/java LL_dup_2.fastq
```

Parameters of FastUniq:
-t q : output sequence format - FASTQ format into TWO output files
-o : first output / -p : second output
-c 1 : types of sequence descriptions for output - new serial numbers assigned by FastUniq (0 : the raw descriptions)

### Mapping

```
mkdir /your-path/PoolSeq_Clec/Mapped/
cd /your-path/PoolSeq_Clec/Mapped/

/your-path/Tools/bwa/bwa mem -t 8 /your-path/PoolSeq_Clec/Ref_Clec/ncbi-genomes-2020-12-15/GCF_000648675.2_Clec_2.1_genomic.fna \
/your-path/PoolSeq_Clec/Trimmed/LL_dup_1.fastq /your-path/PoolSeq_Clec/Trimmed/LL_dup_2.fastq > LL_mapped.sam
```

Now we want to control mapping quality, while converting file from sam to bam, and keeping only reads whose mapping has a probability > 99% to be correct ("-q 20").

Two options: 
- keep all reads (to know the percentage mapping using flagstat)
- separate mapped reads ("-F 4") from unmapped ("-f 4") reads

```
/your-path/Tools/samtools-1.10/bin/samtools view -bS -F 4 -q 20 LL_mapped.sam > LL_mapped.bam
/your-path/Tools/samtools-1.10/bin/samtools view -bS -f 4 -q 20 LL_mapped.sam > LL_unmapped.bam

/your-path/Tools/samtools-1.10/bin/samtools sort LL_mapped.bam -o LL_mapped_sorted.bam

/your-path/Tools/samtools-1.9/bin/samtools flagstat LL_mapped_sorted.bam > flagstat_LL.txt
```

### Analysing coverage 

We can compute coverage per base (i.e. number of reads mapping at a position of the reference genome)

```
/your-path/Tools/Tools/bedtools2/bin/bedtools genomecov -ibam SF_mapped_sorted.bam -d > SF_cov.txt
```

We choose to exclude of our analysis coverage over 50 bp, which corresponds to >95% quantile for all populations, in order to avoid bias due to very high coverage. Coverages computed here were also used for copy number variation analysis.



## Genetic differenciation of populations

The goal was to understand what differenciates the four *Cimex lectularius* PoolSeq samples: London Lab, London Field, German Lab, Sweden Field. 
Our hypothesis was that we could be able to find candidate loci correlated with their insecticide resistance phenotypes - resistant for Field strains and susceptible for Lab strains. 

For following analysis, we excluded the scaffold "NC_030043.1", which corresponds to the mitochondrial genome. Indeed, for one copy of the nuclear genome, there are several copies of the nuclear genome. Furthermore, the mitochondrial genome does not evolve like the nuclear genome (not the same mutation rate, no recombination, maternal transmission). 

### Install tools

Here are the tools and versions used: 
- PoPoolation 2 v1201
- R/packages poolfstat v2.0.0; pcadapt v4.3.3;

They will be store in /your-path/Tools.

### Detecting SNPs

A SNP (Single Nucleotide Polymorphism) is a ponctual mutation that may be associated with candidate regions for resistance.

First, we converted our files in format that can be used by PoPoolation:

``` 
/your-path/Tools/samtools-1.9/bin/samtools index LL_mapped_sorted.bam # do this for all strains

/your-path/samtools-1.9/bin/samtools mpileup -B /your-path/PoolSeq_Clec/Mapped/LL_mapped_sorted.bam /your-path/PoolSeq_Clec/Mapped/LF_mapped_sorted.bam /your-path/PoolSeq_Clec/Mapped/GL_mapped_sorted.bam /your-path/PoolSeq_Clec/Mapped/SF_mapped_sorted.bam > all.mpileup

perl /your-path/Tools/popoolation2_1201/mpileup2sync.pl --fastq-type sanger --min-qual 20 --input all.mpileup --output all.sync
``` 

The output file 'all.sync' can be used with R/poolfstat to compute SNPs:  

``` 
library(poolfstat)

# wc -l all.sync = 509,820,671
# Estimate time - 0.36s/Mi lines processed -> 510*0.36 ~ 184 min soit 3h
# Default parameters - start at 13h45, end at ~17h (189.86 min)

pooldata = popsync2pooldata(sync.file = "all.sync", 
                            poolsizes = rep(30,4), 
                            poolnames = c("GL","LF","SF","LL"), 
                            min.rc = 1,
                            min.cov.per.pool = -1,
                            max.cov.per.pool = 1e+06, 
                            min.maf = 0.01, 
                            noindel = TRUE,
                            nlines.per.readblock = 1e+06, 
                            nthreads = 1)
# 509.8207 millions lines processed in 189.86  min.;  10 139 943 SNPs found for 4 Pools

# Do subset of pooldata :
pooldata_sub <- pooldata.subset(pooldata, pool.index = c(1,2,3,4), 
                                min.cov.per.pool = 10, 
                                max.cov.per.pool = 50, 
                                min.maf = 0.05)
# pool.index = indexes of the pools (at least two), that should be selected to create the new pooldata objec
# min.maf correspond to minimal allowed Minor Allele Frequency
# max.cov of 50 : corresponds to coverage >q95
# 5 574 190 SNPs for 4 pools
``` 

Perform PCA on SNPs:

```
# Build pcadapt matrix

ref_GL <- pooldata_sub@refallele.readcount[,1]
ref_LF <- pooldata_sub@refallele.readcount[,2]
ref_SF <- pooldata_sub@refallele.readcount[,3]
ref_LL <- pooldata_sub@refallele.readcount[,4]

alt_GL <- pooldata_sub@readcoverage[,1] - pooldata_sub@refallele.readcount[,1]
alt_LF <- pooldata_sub@readcoverage[,2] - pooldata_sub@refallele.readcount[,2]
alt_SF <- pooldata_sub@readcoverage[,3] - pooldata_sub@refallele.readcount[,3]
alt_LL <- pooldata_sub@readcoverage[,4] - pooldata_sub@refallele.readcount[,4]

fq_GL <- ref_GL/pooldata_sub@readcoverage[,1]
fq_LF <- ref_LF/pooldata_sub@readcoverage[,2]
fq_SF <- ref_SF/pooldata_sub@readcoverage[,3]
fq_LL <- ref_LL/pooldata_sub@readcoverage[,4]
  
pooldata_sub@snp.info[,1] <- substring(pooldata_sub@snp.info[,1],1,12)
SNP <-paste(pooldata_sub@snp.info[,1],pooldata_sub@snp.info[,2] ,sep="_")

mat_biall_poolfstat=matrix(nrow=4,ncol=5722762)
colnames(mat_biall_poolfstat) <- SNP
rownames(mat_biall_poolfstat)=c("GL","LF","SF","LL")

mat_biall_poolfstat[1,]=fq_GL
mat_biall_poolfstat[2,]=fq_LF
mat_biall_poolfstat[3,]=fq_SF
mat_biall_poolfstat[4,]=fq_LL

mat_biall_poolfstat[1:4,1:10]

#install.packages("pcadapt")
library(pcadapt)

mat_biall_pca <- mat_biall_poolfstat[1:4,]
pca<-read.pcadapt(mat_biall_pca, type = "pool")
res <- pcadapt(pca)
summary(res)

# PC scores
res$scores[,1] # Premier axe
res$scores[,2] # Second axe

hist(res$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
# Confirms that most of the p-values follow an uniform distribution. The excess of small p-values indicates the presence of outliers.

poplist.names <- c("German Lab", "London Field","Sweden Field","London Lab")
plot(res, option = "scores", i = 1, j = 2, pop = poplist.names) # Pour voir PC1 vs PC2
plot(res, option = "manhattan")
```

### Compute FST


For each SNP, we can compute SNP-specific pairwise FST for each comparisons between strains (GL_vs_LF, GL_vs_SF, GL_vs_LL, LF_vs_SF, LF_vs_LL, SF_vs_LL), thanks to the option "output.snp.values = TRUE": 
``` 
PairwiseFST_all = na.omit(computePairwiseFSTmatrix(pooldata_sub, method = "Anova",
						   min.cov.per.pool = -1, 
                                                   max.cov.per.pool = 1e+06,
                                                   min.maf = -1, 
                                                   output.snp.values = TRUE))
``` 



## Evidence of selection

We performed a contrast analysis to identify SNPs associated with populations ecotypes. This trait (populations' ecotype) being binary, we can use C2 statistic (Olazcuaga et al., 2019) to identify those SNPs, rather than parametric models used to estimates Bayes' Factor (BF).

### Installing tools

BayPass :
``` 
/your-path/Tools/BayPass/baypass_2.2/sources/g_baypass
```

### Identifying genetic markers under selection

1st step // Input data for BayPass:
We can convert Poolfstat SNPs data into BayPass input format.

```
pooldata2genobaypass(pooldata_sub,
                     writing.dir = getwd(), # directory where to create the files
                     prefix = "poolfstatdata_220321", # prefix used for output file names
                     subsamplesize = -1) # all SNPs are considered in the output
                     # subsamplingmethod = "thinning")
```

This file contains reads number for each of n SNPs markers of each n sampled populations. It is a matrix with n lines (corresponding to n SNPs) and 2 x n columns (twice the populations number).

Extract:
```
poolfstatdata_220321.genobaypass
3 7 0 7 0 3 0 10
8 2 4 2 4 0 6 5
```
This corresponds to 2 markers and 4 populations: 1st SNP -> 3 copies of 1st allele in 1st population, 7 copies of 2nd allele in 1st population,... 0 copies of 1st allele in 2nd population, 3 copies of 2nd allele in 2nd population,...

We need two other files to make our analysis:

A file containing haploid sizes for pool-seq data:
```
poolfstatdata_110221.poolsize
30 30 30 30
```

And a contrast file, to do analysis in association with binary traits. It allows calculate contrast of standardised allele frequencies between 2 population groups. The group membership for each population is '1' for the first group, '-1' for the alternative group, and '0' if excluded from the contrast analysis.

For populations in this order: German Lab (GL), London Field (LF), Sweden Field (SF), London Lab (LL) ; we want to compute a contrast between Lab and Field populations:
```
poolfstatdata_110221.ecotype
1 -1 -1 1 
```

Since our data are very large (10 139 943 SNPs), we first divided the input file to have 10 000 SNPs/file, i.e. 1014 files:
```
cd /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/
mkdir split_input
split -dl 10000 poolfstatdata_110221.geno --additional-suffix=.geno /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/split_input/poolfstatdata_110221_
```

For each input, a corresponding job file is created:
```
mkdir split_output

python
>>> 
import os 
import re
for filename in os.listdir("/beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/split_input/"):
        if filename.endswith(".geno"): 
          nb=re.sub('.*110221_','',filename)
          nb=re.sub('.geno','',nb)
          w = open("/beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/split_input/Cimex_BayPassjob_"+nb+".sh",'w')
          w.write("#!/bin/bash\n")
          w.write("#SBATCH --cpus-per-task=1\n")
          w.write("#SBATCH --mem 10G\n")
          w.write("#SBATCH --exclude=pbil-deb[27]\n")
          w.write("#SBATCH -t 2:00:00\n")
          w.write("#SBATCH -e /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/split_output/baypass_"+nb+".error\n")
          w.write("#SBATCH -o /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/split_output/baypass_"+nb+".out\n")
          w.write("#SBATCH -J Cimex_BaypassJob_"+nb+"\n")
          w.write("date;hostname;pwd\n")
          w.write("baypass=/beegfs/data/soft/BayPass/baypass_2.2/sources/g_baypass\n")
          w.write("$baypass -gfile /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/split_input/poolfstatdata_110221_"+nb+".geno -poolsizefile /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/baypass_input/poolfstatdata_110221.poolsize -d0yij 6 -contrastfile /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/baypass_input/poolfstatdata_110221.ecotype -efile /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/baypass_input/poolfstatdata_110221.ecotype -outprefix poolfstatdata_110221_"+nb+"")
          w.write("date")
```

The initial delta (δ) of the distribution of the yij proposal (-d0yij parameter) is generally 1/5 of the size of the smallest pool: 30/5=6. We have followed BayPass pipeline for pool-seq data.
Estimate time running for BayPass on each sub-file is approximately 45 minutes. We can launch all jobs simultaneously using:
```
cd /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/split_input/
for file in ls *sh
do
    sbatch $file
done
```

Before concatenating the outputs, the first line (column names) of each file must be remove:
```
sed -i '1d' *constrast.out
sed -i '1d' *_summary_betai_reg.out
```

Concatenate the outputs:
```
ls -v poolfstatdata_110221_*_summary_contrast.out | xargs cat > poolfstatdata_110221_all_summary_contrast.out 
ls -v poolfstatdata_110221_*_summary_betai_reg.out | xargs cat > poolfstatdata_110221_all_summary_betai_reg.out 
# To make sure that the files are in "numerical" order
```

Check that the two files we want to paste have the same number of lines:
```
wc -l /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/baypass_input/poolfstatdata_110221.snpdet # 10 139 943
wc -l poolfstatdata_110221_all_summary_contrast.out # 10 139 943
wc -l poolfstatdata_110221_all_summary_betai_reg.out # 10 139 943
```

Add column names:
```
sed  -i '1i SCAFFOLD POSITION REF ALT' /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/baypass_input/poolfstatdata_110221.snpdet
sed  -i '1i CONTRAST MRK M_C2 SD_C2 C2_std log10(1/pval)' /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/split_output/poolfstatdata_110221_all_summary_contrast.out
sed  -i '1i COVARIABLE MRK M_Pearson SD_Pearson BF(dB) Beta_is SD_Beta_is eBPis' /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/split_output/poolfstatdata_110221_all_summary_betai_reg.out
```

Paste tables:
```
paste /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/baypass_input/poolfstatdata_110221.snpdet poolfstatdata_110221_all_summary_contrast.out > poolfstatdata_110221_all_summary_contrast_snpdet.out
paste /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/baypass_input/poolfstatdata_110221.snpdet poolfstatdata_110221_all_summary_betai_reg.out > poolfstatdata_110221_all_summary_betai_reg_snpdet.out
```

Open R on the cluster:
```
/beegfs/data/soft/R-3.5.2/bin/R
contrast=read.table("poolfstatdata_110221_all_summary_contrast_snpdet.out")
contrast <- contrast[-1, ]
contrast <- contrast[,-c(5,6)] # Remove column CONTRAST (only 1) and MRK (1:10000 for each subset file)
colnames(contrast) <- c("SCAFFOLD","POSITION","REF","ALT","M_C2","SD_C2",
"C2_std","LOG") # LOG=log10(1/pval)
dim(contrast) # [1] 10139943        8
contrast <- contrast[!(contrast$SCAFFOLD=="NC_030043.1"),] # Remove mitochondrial
dim(contrast) #[1] 10139835        8
contrast$LOG <- as.numeric(as.character(contrast$LOG))

betai=read.table("poolfstatdata_110221_all_summary_betai_reg_snpdet.out")
betai <- betai[-1, ]
betai <- betai[,-c(5,6)] # Remove column CONTRAST (only 1) and MRK (1:10000 for each subset file)
colnames(betai) <- c("SCAFFOLD","POSITION","REF","ALT","M_Pearson",
"SD_Pearson","BF.dB","Beta_is","SD_Beta_is","eBPis") # BF.dB = BF(dB)
dim(betai) # [1] 10139943        8
betai <- betai[!(test3$SCAFFOLD=="NC_030043.1"),] # Remove mitochondrial
dim(betai) #[1] 10139835       10
betai$BF.dB <- as.numeric(as.character(betai$BF.dB))

summary(contrast)
         M_C2                 SD_C2                 C2_std        
 0.68335906:       5   0.87546325:       5   0.00000000:     531  
 0.68453837:       5   0.93783286:       5   0.00000001:     400  
 0.70137523:       5   0.94074307:       5   0.00000002:     279  
 0.70264553:       5   0.94538022:       5   0.00000003:     220  
 0.71227663:       5   1.00928620:       5   0.00000004:     189  
 0.71717973:       5   1.01447543:       5   0.00000005:     161  
 (Other)   :10139913   (Other)   :10139913   (Other)   :10138163  
      LOG        
 Min.   :0.0000  
 1st Qu.:0.1398  
 Median :0.3639  
 Mean   :0.4971  
 3rd Qu.:0.7124  
 Max.   :8.4857  

summary(betai) 
       M_Pearson             SD_Pearson           BF.dB        
 -0.00120966:       5   0.54029289:       9   Min.   :-16.755  
 -0.00761172:       5   0.51934770:       8   1st Qu.: -7.034  
 -0.01160689:       5   0.52120184:       8   Median : -5.397  
 -0.01261898:       5   0.52303979:       8   Mean   : -5.139  
 -0.02761936:       5   0.52951295:       8   3rd Qu.: -3.858  
 -0.04455849:       5   0.53196095:       8   Max.   : 59.787  
 (Other)    :10139805   (Other)   :10139786                    
        Beta_is              SD_Beta_is              eBPis         
 -0.00064555:      18   0.01555728:      17   0.00321555:       7  
 0.00277646 :      18   0.01961458:      17   0.03123584:       7  
 0.00074518 :      17   0.01459807:      16   0.04042228:       7  
 -0.00475267:      16   0.01537673:      16   0.07455867:       7  
 0.00052373 :      16   0.01680086:      16   0.08742999:       7  
 0.00082491 :      16   0.01717066:      16   0.00012522:       6  
 (Other)    :10139734   (Other)   :10139737   (Other)   :10139794  
```

We are looking for markers with C2 value significantly different from 0 (low p-value), which means that those markers are associated with the population ecotype (here, field vs lab strains). Since Bayes Factor (BF) measures the likelihood of a model under selection, we also track high BF.
The resulting C2 contrasts (and BF) might then be plotted (and compared) as follows:
```
#check the behavior of the p-values associated to the C2
hist(10**(-1*test$LOG),freq=F,breaks=100)
abline(h=1, col="red")

plot(betai$BF.dB, test$LOG, xlab="BF",ylab="C2 p-value (-log10 scale)", pch=20)
abline(h=3,lty=2,col="blue") #0.001 p--value theshold
abline(v=20,lty=2,col="red") #BF threshold for decisive evidence (according to Jeffreys' rule)
```

<img src="plot_poolfstatdata_110221_hist.png" style="background:none; border:none; box-shadow:none;">

The histogram representing the associated p-values is not regular.

How many selected points are in the top right corner of the plot?
```
all <- merge(x = contrast[ , c("SCAFFOLD","POSITION","LOG")], y = betai[ , c("SCAFFOLD","POSITION","BF.dB")], by=c("SCAFFOLD","POSITION"),)
write.table(all, file="/beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/baypass_110221_results.txt")

outliers <- all[(all$BF.dB>20 & all$LOG>3),]
dim(outliers) # [1] 261   4
write.table(outliers, file="/beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/baypass_110221_outliers.txt")
```
Here we have selected 261 SNPs that might be subjected to selection according to our BayPass analysis between field and laboratory populations.

We now use the sub-poolfstat data, with more stringent parameters: min.cov.per.pool = 10 (previously = -1), max.cov.per.pool = 200 (previously = 1e+06), min.maf = 0.05 (previously = 0.01).
Note: we tried to subset the result table, but LOG values were different from when running BayPass with a subset input SNPs file.

```
mkdir split_input_sub
split -dl 10000 poolfstatdata_sub_110221.geno --additional-suffix=.geno /beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/split_input_sub/poolfstatdata_sub_110221_
```

As before, we have created jobs for each input files, and launching them all simultaneously. The first line (column name) of each ouput file with 10000 lines must be removed before concatenating the outputs. Output files are then pasted with snps_det fail, containing all informations about SNP positions and scaffold name. Files are then processed with R as above (on data not undersampled).

The resulting C2 contrasts (and BF) might then be plotted (and compared) as follows:
```
#check the behavior of the p-values associated to the C2
hist(10**(-1*contrast$LOG),freq=F,breaks=100)
abline(h=1, col="red")

plot(betai$BF.dB, contrast$LOG, xlab="BF",ylab="C2 p-value (-log10 scale)", pch=20)
abline(h=3,lty=2,col="blue") #0.001 p--value theshold
abline(v=20,lty=2,col="red") #BF threshold for decisive evidence (according to Jeffreys' rule)
```
<img src="plot_poolfstatdata_sub_110221.png" style="background:none; border:none; box-shadow:none;">

We can observe an histogramm with anti-conservative p-values, which seems much more reliable.

How many selected points are in the top right corner of the plot?
```
all <- merge(x = contrast[ , c("SCAFFOLD","POSITION","LOG")], y = betai[ , c("SCAFFOLD","POSITION","BF.dB")], by=c("SCAFFOLD","POSITION"),)
write.table(all, file="/beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/baypass_sub_110221_results.txt")

outliers <- all[(all$BF.dB>20 & all$LOG>3),]
dim(outliers) # 168
write.table(outliers, file="/beegfs/data/chaberkorn/PoolSeq_Clec/BayPass/baypass_sub_110221_outliers.txt")
```
Here we have selected 168 SNPs that might be subjected to selection according to our BayPass analysis between field and laboratory populations.


### Mapping areas under selection

We now need to map those outliers on linkage group (LG) and without LG, to see if their positions coincide with the presence of high FSTs between resistant and susceptible strains. We only work on sub_outliers_wo_LG as p-values distribution seem to show less bias.
First, we do not use LG :

```
sub_outliers <- read.table(file="~/Desktop/CloudStation/THESE/WholeGenome PoolSeq/BayPass/baypass_sub_110221_outliers.txt", 
                           sep=" ", header=T)
colnames(sub_outliers) <- c("scaffold","position","LOG","BF.dB")
sub_outliers[,1] <- substring(sub_outliers[,1],1,12)

load("/Users/chloe/Documents/Cluster/Envir_FST_Poolfstat.RData")

# Since we are working on outliers detected on pooldata_sub, we shoud use these files :
dim(PairwiseFST_all$PairwiseSnpFST)
dim(pooldata_sub@snp.info)
# Verify same dim()

data_FST <- cbind(pooldata_sub@snp.info,PairwiseFST_all$PairwiseSnpFST)
data_FST <- as.data.frame(data_FST)
colnames(data_FST) <-c("scaffold","position","ref","alt","GL_vs_LF","GL_vs_SF",
                       "GL_vs_LL","LF_vs_SF","LF_vs_LL","SF_vs_LL")
# Clean our table :
data_FST[,1] <- substring(data_FST[,1],1,12) # 5 722 762
data_FST <- data_FST[!(data_FST$scaffold=="NC_030043.1"),] # 5 722 749
data_FST <- data_FST[,-c(3,4)] # enlever col ref et alt

summary(data_FST) # warning : all variables are characters
data_FST$position <- as.numeric(data_FST$position)
data_FST$GL_vs_LF <- as.numeric(data_FST$GL_vs_LF)
data_FST$GL_vs_SF <- as.numeric(data_FST$GL_vs_SF)
data_FST$GL_vs_LL <- as.numeric(data_FST$GL_vs_LL)
data_FST$LF_vs_SF <- as.numeric(data_FST$LF_vs_SF)
data_FST$LF_vs_LL <- as.numeric(data_FST$LF_vs_LL)
data_FST$SF_vs_LL <- as.numeric(data_FST$SF_vs_LL)

data_FST$Colour="black"
data_FST$Colour[(data_FST$GL_vs_LF > quantile(data_FST$GL_vs_LF, 0.99, na.rm=T) 
                   & data_FST$GL_vs_SF > quantile(data_FST$GL_vs_SF, 0.99, na.rm=T)
                   & data_FST$LF_vs_LL > quantile(data_FST$LF_vs_LL, 0.99, na.rm=T)
                   & data_FST$SF_vs_LL > quantile(data_FST$SF_vs_LL, 0.99, na.rm=T))]="red"

nb_FST <- data_FST[data_FST$Colour == "red",] # 231 FST R≠S q99

# Create a column to know which lines are BayPass outliers :
sub_outliers$out <- "yes" 
sub_outliers_wo_LG <- merge(x = data_FST[ , c("scaffold","GL_vs_LF","GL_vs_SF","GL_vs_LL","LF_vs_SF","LF_vs_LL","SF_vs_LL","position","Colour")], 
                         y = sub_outliers[ , c("scaffold","position","out")], by=c("scaffold","position"),all=T)

library(dplyr)
sub_outliers_wo_LG$X <- 1:length(sub_outliers_wo_LG$scaffold)
sub_outliers_wo_LG <- sub_outliers_wo_LG %>% arrange(position)
sub_outliers_wo_LG <- sub_outliers_wo_LG %>% arrange(scaffold)
sub_outliers_wo_LG <- sub_outliers_wo_LG %>% arrange(Colour)

View(sub_outliers_wo_LG)

nb_FST_out <- subset(sub_outliers_wo_LG , sub_outliers_wo_LG$Colour == "red" & sub_outliers_wo_LG$out == "yes")
# 16 SNPs combine both criteria

library(scales)
cols <- sub_outliers_wo_LG$Colour
plot(x=sub_outliers_wo_LG$X,y=sub_outliers_wo_LG$LF_vs_LL, ylab="", xlab="", pch=20, main="FST LF_vs_LL", col=alpha(cols, 0.4))
points(x=sub_outliers_wo_LG$X[sub_outliers_wo_LG$out == "yes"], 
       y=sub_outliers_wo_LG$LF_vs_LL[sub_outliers_wo_LG$out == "yes"],
       pch = 20, col =alpha("blue", 0.4))
```

<img src="plot_FST_baypass_in_blue.png" width="792" height="546" style="background:none; border:none; box-shadow:none;">


