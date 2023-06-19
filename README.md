![Image 1](bedbugs.png)

You can learn more about the context of this repository in the article from which it stems: https://onlinelibrary.wiley.com/doi/full/10.1111/eva.13550

Contact: chloe.haberkorn@univ-lyon1.fr / chloehbk@gmail.com

## Table of Contents

- **[Installing tools](#Installing-tools)**

- **[Pool-seq data processing](#Pool-seq-data-processing)**
	- [Getting the data](#Getting-the-data)
	- [Trimming](#Trimming)
	- [Removing duplicates](#Removing-duplicates)
	- [Mapping](#Mapping)
	- [Analysing coverage](#Analysing-coverage)

- **[Overall SNPs analyzes](#Overall-SNPs-analyzes)**
	- [Detecting SNPs](#Detecting-SNPs)
	- [Performing PCA](#Performing-PCA)
	- [Computing FST](#Computing-FST)

 - **[Selecting candidate SNPs](#Selecting-candidate-SNPs)**
	- [Differentiated FST](#Differentiated-FST)
	- [Selection with contrast between phenotypes](#Selection-with-contrast-between-phenotypes)
	- [Subsetting data on Minor Allele Frequency](#Subsetting-data-on-Minor-Allele-Frequency)
	- [Alternative alleles](#Alternative-alleles)
	- [Combining conditions](#Combining-conditions)
	- [Synonymous or not](#Synonymous-or-not)

 - **[Copy Number Variation](#Copy-Number-Variation)**
 	- [Computing average coverage by gene](#Computing-average-coverage-by-gene)
	- [Selecting amplified genes](#Selecting-amplified-genes)



## Installing tools

Here are the tools and versions used (on a computing cluster): 
- FastQC 
- Trimmomatic v0.39
- FastUniq v1.1
- BWA v0.74
- Samtools v1.9
- Bedtools v2.29.1
- PoPoolation 2 v1201
- BayPass v2.3
- R v3.5.2

They will be store in /your-path/Tools.
We also used R on a computer with packages poolfstat v2.0.0, pcadapt v4.3.3, VariantAnnotation v1.34.0, and GenomicFeatures v1.40.1.


## Pool-seq data processing

The goal is first to map *Cimex lectularius* PoolSeq samples (London Lab, London Field, German Lab and Sweden Field - pools of 30 individuals) on reference genome.

### Getting the data

Raw sequences (fastq.gz files) are available on SRA under the accession number PRJNA826750:
- London Lab (LL): SRX15014458
- London Field (LF): SRX15014459
- German Lab (GL): SRX15014460
- Sweden Field (SF): SRX15014461

```
mkdir /your-path/PoolSeq_Clec/Raw_Clec
```
Please keep in mind that the following step should be performed for each strains, but will be presented for only one for the sake of brevity.

We used the recent reference genome and annotation, avalaible here: https://www.ncbi.nlm.nih.gov/assembly/GCF_000648675.2
```
mkdir /your-path/PoolSeq_Clec/Ref_Clec
```

### Trimming

To use Trimmomatic, we also had to download the associated adapters. To know which adapters should be downloaded, one can look at overrepresented sequences in fastq file: in our case, "TruSeq Adapter" > https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE.fa

We kept Trimmomatic default parameters, except for:
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

We first compute coverage per base (i.e. number of reads mapping at a position of the reference genome):

```
/your-path/Tools/Tools/bedtools2/bin/bedtools genomecov -ibam SF_mapped_sorted.bam -d > SF_basecov.txt
```

We choose to exclude of our analysis coverage over 50 bp, which corresponds to >95% quantile for all populations, in order to avoid bias due to very high coverage. We were then able to compute the coverage (X) for the 4 strains:

```
awk '{ total += $3 } END { print total/NR }' SF_basecov.txt 
```
Coverage LL = 25.0468
Coverage LF = 32.0371
Coverage GL = 39.4606
Coverage SF = 25.3821

## Overall SNPs analyzes

The goal was to understand what differenciates the four *Cimex lectularius* PoolSeq samples: London Lab, London Field, German Lab, Sweden Field. 
Our hypothesis was that we could be able to find candidate loci correlated with their insecticide resistance phenotypes - resistant for Field strains and susceptible for Lab strains. 

For following analysis, we excluded the scaffold "NC_030043.1", which corresponds to the mitochondrial genome. Indeed, for one copy of the nuclear genome, there are several copies of the nuclear genome. Furthermore, the mitochondrial genome does not evolve like the nuclear genome (not the same mutation rate, no recombination, maternal transmission). 

### Detecting SNPs

A SNP (Single Nucleotide Polymorphism) is a ponctual mutation that may be associated with candidate regions for resistance.

First, we converted our files in format that can be used by PoPoolation:

``` 
/your-path/Tools/samtools-1.9/bin/samtools index LL_mapped_sorted.bam # do this for all strains

/your-path/samtools-1.9/bin/samtools mpileup -B /your-path/PoolSeq_Clec/Mapped/GL_mapped_sorted.bam /your-path/PoolSeq_Clec/Mapped/LF_mapped_sorted.bam /your-path/PoolSeq_Clec/Mapped/SF_mapped_sorted.bam /your-path/PoolSeq_Clec/Mapped/LL_mapped_sorted.bam > all.mpileup

perl /your-path/Tools/popoolation2_1201/mpileup2sync.pl --fastq-type sanger --min-qual 20 --input all.mpileup --output all.sync
``` 

SNPs can be detected on the output file 'all.sync' using R/poolfstat:  

``` 
library(poolfstat)

# wc -l all.sync = 509,820,671 > 1 ligne par position du génome
# Poolsizes correspond to haploid size (30*2)
# Here, creating a job on cluster with 8 threads

pooldata = popsync2pooldata(sync.file = "all.sync", 
                            poolsizes = c(60,60,60,56), 
                            poolnames = c("GL","LF","SF","LL"),
                            min.rc = 1,
                            min.cov.per.pool = -1,
                            max.cov.per.pool = 1e+06, 
                            min.maf = 0.01, 
                            noindel = TRUE,
                            nlines.per.readblock = 1e+06, 
                            nthreads = 8)
# 509.8207 millions lines processed in 41.8  min.;  10,139,943 SNPs found SNPs found for 4 Pools

# Do subset of pooldata :
pooldata_sub <- pooldata.subset(pooldata, pool.index = c(1,2,3,4), 
                                min.cov.per.pool = 10, 
                                max.cov.per.pool = 50, 
                                min.maf = 0.01)
# pool.index = indexes of the pools (at least two), that should be selected to create the new pooldata object
# min.maf correspond to minimal allowed Minor Allele Frequency
# max.cov of 50 : corresponds to coverage >q95
# Data consists of 8,034,712 SNPs (8,03 M) SNPs for 4 Pools
``` 

### Performing PCA

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

mat_biall_poolfstat=matrix(nrow=4,ncol=length(pooldata_sub@snp.info))
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
res$singular.values
# singular.values is a vector containing the K ordered square root of the proportion of variance explained by each PC.
PC1=res$singular.values[1]^2
PC2=res$singular.values[2]^2

poplist.names <- c("German Lab", "London Field","Sweden Field","London Lab")
plot(res, option = "scores", i = 1, j = 2, pop = poplist.names) # PC1 vs PC2
```

### Computing FST

For each SNP, we can compute SNP-specific pairwise FST for each comparisons between strains (GL_vs_LF, GL_vs_SF, GL_vs_LL, LF_vs_SF, LF_vs_LL, SF_vs_LL), thanks to the option "output.snp.values = TRUE": 
``` 
PairwiseFST_all = na.omit(compute.pairwiseFST(pooldata_sub, method = "Anova",
					      min.cov.per.pool = -1, 
                                              max.cov.per.pool = 1e+06,
                                              min.maf = -1, 
                                              output.snp.values = TRUE))
# Warning : here, min.maf by pairs, not for the 4 strains together !

head(PairwiseFST_all@values) # to build Table 2
``` 


## Selecting candidate SNPs

### Differentiated FST

We first merged informations together for following analysis:
``` 
# append  allele frequencies
ps=pooldata_sub2@refallele.readcount/pooldata_sub2@readcoverage

# construct a table containing all useful info
data=data.frame(pooldata_sub2@snp.info, 
                pooldata_sub2@refallele.readcount, 
                pooldata_sub2@readcoverage, 
                ps, 
                PairwiseFST_all@PairwiseSnpFST)

names(data)=c("contig", "position", "ref", "alt", 
              paste0(pooldata_sub2@poolnames, "_ref"),
              paste0(pooldata_sub2@poolnames, "_tot"), 
              paste0(pooldata_sub2@poolnames, "_p"),
              "GL_vs_LF", "GL_vs_SF", "GL_vs_LL","LF_vs_SF" ,"LF_vs_LL","SF_vs_LL")

data <- data[!(data$contig=="NC_030043.1"),] # Remove mitochondrial
``` 

### Selection with contrast between phenotypes

We performed a contrast analysis to identify SNPs associated with populations ecotypes. This trait (populations' ecotype) being binary, we can use C2 statistic (Olazcuaga et al., 2019) to identify those SNPs, rather than parametric models used to estimates Bayes' Factor (BF).

First, we converted Poolfstat SNPs data into BayPass input format:
```
pooldata2genobaypass(pooldata_sub,
                     writing.dir = getwd(), # directory where to create the files
                     prefix = "poolfstatdata_080422", # prefix used for output file names
                     subsamplesize = -1) # all SNPs are considered in the output
                     # subsamplingmethod = "thinning")
```

This created 3 files. The first one contains reads number for each of n SNPs markers of each p sampled populations: it's a matrix with n lines (corresponding to n SNPs) and 2 x p columns (twice the populations number).

Extract:
```
poolfstatdata_080422.genobaypass
3 7 0 7 0 3 0 10
8 2 4 2 4 0 6 5
```
This corresponds to 2 markers and 4 populations

1st line : 1st SNP (marker) -> 3 copies of 1st allele in 1st population ; 7 copies of 2nd allele in 1st population ; 0 copies of 1st allele in 2nd population ; 7 copies of 2nd allele in 2nd population ; ...

We need two other files to make our analysis:

A file containing haploid sizes for pool-seq data:
```
poolfstatdata_080422.poolsize
60 60 60 56
```

And a contrast file, to do analysis in association with binary traits. It allowed to compute contrast of standardized allele frequencies between 2 population groups. The group membership for each population is '1' for the first group, '-1' for the alternative group, and '0' if excluded from the contrast analysis.

For populations in this order: German Lab (GL), London Field (LF), Sweden Field (SF), London Lab (LL). We want to compute a contrast between Lab and Field populations:
```
poolfstatdata_080422.ecotype
1 -1 -1 1 
```

Running Baypass:
```
mkdir /your-path/PoolSeq_Clec/BayPass/
cd /your-path/PoolSeq_Clec/BayPass/

baypass=/your-path/Tools/BayPass/baypass_2.3/sources/g_baypass
$baypass -gfile /your-path/PoolSeq_Clec/BayPass/poolfstatdata_080422.genobaypass \
-poolsizefile /your-path/PoolSeq_Clec/BayPass/poolfstatdata_080422.poolsize -d0yij 6 \
-contrastfile /your-path/PoolSeq_Clec/BayPass/baypass_input/poolfstatdata_080422.ecotype \
-efile /your-path/PoolSeq_Clec/BayPass/poolfstatdata_080422.ecotype -outprefix poolfstatdata_080422 -nthreads 16
```

The initial delta (δ) of the distribution of the yij proposal (-d0yij parameter) is generally 1/5 of the size of the smallest pool: 30/5=6. We followed BayPass pipeline for pool-seq data.

We were looking for markers with C2 value significantly different from 0 (low p-value), which means that those markers are associated with the population ecotype (here, field vs lab strains). Since Bayes Factor (BF) measures the likelihood of a model under selection, we also tracked high BF.

We then merged (with R for example) two of output files together: poolfstatdata_080422_summary_contrast.out.out and poolfstatdata_080422_summary_betai_reg.out, in /your-path/PoolSeq_Clec/BayPass/baypass_080422_results.txt")

We built a table combining all informations:
``` 
baypass=read.table("/your-path/baypass_080422_results.txt", sep=" ")
colnames(baypass) <- c("contig","position","C2_std","LOG_C2","BF.dB","XtXst","LOG_xtx")

data <- merge(data, baypass, by.x=c("contig","position"), by.y=c("SCAFFOLD","POSITION"), all.y=F) # add baypass to SNP/FST informations
``` 

Since reference allele is "chosen arbitrarily in each pool" by poolfstat, we need to correct it. Indeed, we want the reference allele to match the reference genome allele, i.e. Harlan strain allele.

#1 Add a column with reference allele related to Harlan, and correct ref/alt columns

```{r}
data_bed <- data[,c(1,2)]
data_bed$position <- data_bed$position-1

library(dplyr)
data_bed = data_bed %>% arrange(position)

write.table(data_bed, file = "~/data_bed.txt",quote=F, row.names=F, col.names=F, sep="\t")
```

Then, use this tool: https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html
```
/beegfs/data/chaberkorn/Tools/bedtools2/bin/bedtools getfasta -fi /beegfs/data/chaberkorn/PoolSeq_Clec/GCF_000648675.2_Clec_2.1_genomic.fna -bed data.bed -fo data.fa.out
```
! Check that there is no number in scientific writing in BED entry ! 

Then, we created a column in our "data" table called "nucleotides" which corresponds to the reference nucleotide extracted from "data.fa.out".
This column must be formated just as our "ref" and "alt" columns:
```{r}
data$nucleotides <- str_to_upper(data$nucleotides) # all nucleotides in upper case

# Check for differences between Poolfstat "ref" allele and reference genome "nucleotides"
wrong_ref <- (data$ref != data$nucleotides) & (data$alt ==  data$nucleotides)
wrong_ref_count <- data[(data$ref != data$nucleotides) & (data$alt ==  data$nucleotides),] 

true_ref <- (data$ref == data$nucleotides) & (data$alt !=  data$nucleotides)

# Check whether there are alleles of ref & alt != alleles of the reference genome (harlan strain)
wrong_all <- (data$ref != data$nucleotides) & (data$alt !=  data$nucleotides)
# yes, indeed!

# Replace the count of ref by the count of alt 
library(data.table)
data <- as.data.table(data)
data <- data[wrong_ref, LL_ref := LL_tot-LL_ref]
data <- data[wrong_ref, GL_ref := GL_tot-GL_ref]
data <- data[wrong_ref, LF_ref := LF_tot-LF_ref]
data <- data[wrong_ref, SF_ref := SF_tot-SF_ref]

# Replace alternative allele by Poolfstat reference allele and Poolfstat reference allele by Harlan reference allele
data <- data[wrong_ref, alt := ref]
data <- data[wrong_ref, ref := nucleotides]
```

#2 Modify reference allele frequencies (p)

```{r}
data$GL_p <- data$GL_ref/data$GL_tot
data$LL_p <- data$LL_ref/data$LL_tot
data$LF_p <- data$LF_ref/data$LF_tot
data$SF_p <- data$SF_ref/data$SF_tot
```

### Subsetting data on Minor Allele Frequency

We used the BayPass software to estimate the nucleotide diversity "π" for each SNP (see [Selection with contrast between phenotypes](#Selection-with-contrast-between-phenotypes) to know how to run BayPass). With π, we were then able to compute "MAF", the minor allele frequency:

``` 
pi=read.table("poolfstatdata_080422_summary_pi_xtx.out",header=T)
pi$MAF <- 0.5-abs(result$M_P-0.5)
summary(pi$M_P) # Min 0.02142 - Max 0.97752 
summary(pi$MAF) # Min 0.02142 - Max 0.5

pi_sub <- pi[pi$MAF>0.2,]
``` 

Subsetting our data on MAF>0.2 gave us 2,918,170 SNPs (without the mitochondrial scaffold).

### Combining conditions

With BayPass results, we needed to compute adjusted p-values for C2 statistic:

```
data$C2_pval <- 1/(10^(data$LOG_C2))

library(stringi)
library(qvalue)

qobj <- qvalue(p = data$C2_pval)
lfdr <- qobj$lfdr
summary(lfdr)
hist(lfdr, breaks=100)

lfdr_signif_0.2 <- lfdr[lfdr<0.2] # 1,052
```

We want to identify SNPs where London Lab carries Harlan reference allele and London Field an alternative allele, with LL_p and LF_p the frequencies of reference allele.
We also used a threshold of BF > 5, and FST empirical p-values between London Lab and London Field < 0.05:

```
data$Colour="black"
data$Colour[(data$BF.dB>5 
             & data$lfdr_C2<0.2
             & data$Fst_LLLF_pval<0.05
             & data$LL_p>data$LF_p)
            ]="red"
summary(as.factor(data$Colour))
```

By doing so, we detected 580 candidates.

### Synonymous or not






