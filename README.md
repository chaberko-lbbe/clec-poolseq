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
	- [Identifying differentiated SNPs (FST)](#Identifying-differentiated-alleles)
	- [SNPs under selection with contrast between phenotypes (BayPass)](#SNPs-under-selection-with-contrast-between-phenotypes)
	- [Selecting alternative alleles](#Selecting-alternative-alleles)

## Pool-seq data processing

The goal is first to map *Cimex lectularius* PoolSeq samples (London Lab, London Field, German Lab and Sweden Field - pools of 30 individuals) on reference genome.
We used the recent reference genome and annotation, avalaible here: https://www.ncbi.nlm.nih.gov/assembly/GCF_000648675.2

We will have to download a few softs.

### Installing tools

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





## Detecting Single Nucleotide Polymorphism

The goal was to understand what differenciates the four *Cimex lectularius* PoolSeq samples: London Lab, London Field, German Lab, Sweden Field. 
Our hypothesis was that we could be able to find candidate loci correlated with their insecticide resistance phenotypes - resistant for Field strains and susceptible for Lab strains. 

For following analysis, we excluded the scaffold "NC_030043.1", which corresponds to the mitochondrial genome. Indeed, for one copy of the nuclear genome, there are several copies of the nuclear genome. Furthermore, the mitochondrial genome does not evolve like the nuclear genome (not the same mutation rate, no recombination, maternal transmission). 

### Installing tools

Here are the tools and versions used: 
- PoPoolation 2 v1201
- R/packages poolfstat v2.0.0; pcadapt v4.3.3;

They will be store in /your-path/Tools.

### Overall SNPs analyzes

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
For each SNP, we can compute SNP-specific pairwise FST for each comparisons between strains (GL_vs_LF, GL_vs_SF, GL_vs_LL, LF_vs_SF, LF_vs_LL, SF_vs_LL), thanks to the option "output.snp.values = TRUE": 
``` 
PairwiseFST_all = na.omit(computePairwiseFSTmatrix(pooldata_sub, method = "Anova",
						   min.cov.per.pool = -1, 
                                                   max.cov.per.pool = 1e+06,
                                                   min.maf = -1, 
                                                   output.snp.values = TRUE))
``` 




## Selecting candidate SNPs


### Installing tools

BayPass :
``` 
/your-path/Tools/BayPass/baypass_2.2/sources/g_baypass
```

### Identifying differentiated SNPs (FST)

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
This corresponds to 2 markers and 4 populations

1st line : 1st SNP (marker) -> 3 copies of 1st allele in 1st population ; 7 copies of 2nd allele in 1st population ; 0 copies of 1st allele in 2nd population ; 7 copies of 2nd allele in 2nd population ; ...

We need two other files to make our analysis:

A file containing haploid sizes for pool-seq data:
```
poolfstatdata_220321.poolsize
30 30 30 30
```

And a contrast file, to do analysis in association with binary traits. It allowed to compute contrast of standardized allele frequencies between 2 population groups. The group membership for each population is '1' for the first group, '-1' for the alternative group, and '0' if excluded from the contrast analysis.

For populations in this order: German Lab (GL), London Field (LF), Sweden Field (SF), London Lab (LL). We want to compute a contrast between Lab and Field populations:
```
poolfstatdata_220321.ecotype
1 -1 -1 1 
```

Running Baypass:
```
mkdir /your-path/PoolSeq_Clec/BayPass/
cd /your-path/PoolSeq_Clec/BayPass/

baypass=/your-path/Tools/BayPass/baypass_2.2/sources/g_baypass
$baypass -gfile /your-path/PoolSeq_Clec/BayPass/poolfstatdata_220321"+nb+".geno -poolsizefile /your-path/PoolSeq_Clec/BayPass/poolfstatdata_220321.poolsize -d0yij 6 -contrastfile /your-path/PoolSeq_Clec/BayPass/baypass_input/poolfstatdata_220321.ecotype -efile /your-path/PoolSeq_Clec/BayPass/poolfstatdata_220321.ecotype -outprefix poolfstatdata_220321
```

The initial delta (δ) of the distribution of the yij proposal (-d0yij parameter) is generally 1/5 of the size of the smallest pool: 30/5=6. We followed BayPass pipeline for pool-seq data.

We were looking for markers with C2 value significantly different from 0 (low p-value), which means that those markers are associated with the population ecotype (here, field vs lab strains). Since Bayes Factor (BF) measures the likelihood of a model under selection, we also tracked high BF.

We then merged two of output files together: poolfstatdata_220321_summary_contrast_snpdet.out and poolfstatdata_220321_summary_betai_reg_snpdet.out
> /your-path/PoolSeq_Clec/BayPass/baypass_220321_results.txt")



### SNPs under selection with contrast between phenotypes (BayPass)

We performed a contrast analysis to identify SNPs associated with populations ecotypes. This trait (populations' ecotype) being binary, we can use C2 statistic (Olazcuaga et al., 2019) to identify those SNPs, rather than parametric models used to estimates Bayes' Factor (BF).


### Selecting alternative alleles




