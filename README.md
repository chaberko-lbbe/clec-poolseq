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

 - **[Detection of structural variants](#Detection-of-structural-variants)**
 	- [Subset files on scaffolds](#Subset-files-on-scaffolds)
	- [Identify distant and everted read-pairs](#Identify-distant-and-everted-read-pairs)
	- [Cluster distant and everted read-pairs within pools](#Cluster-distant-and-everted-read-pairs-within-pools)
	- [Calculate relative read depth differences](#Calculate-relative-read-depth-differences)
	- [Identify distant and everted read-pairs](#Identify-distant-and-everted-read-pairs)
	- [Append additional allele frequency information](#Append-additional-allele-frequency-information)


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

## Detection of structural variants

We aim to detect two types of structural variants (SVs), duplications and inversions. London Field population was compared to the susceptible reference genome, and to the susceptible London Lab population, in order to detect SVs that could have been selected and underlie the resistance phenotype.
The presence of SVs was inferred based on abnormal read pair orientation and/or distance (insert size), and read depth variation. We used poolCNVcomp (North et al., 2020; Schrider et al., 2013; Schrider et al., 2016), a series of scripts designed for the detection of tandem duplications in pool-seq data, and adapted them to detect inverted duplications and simple inversions as well. 
Associated python scripts can also be downloaded here: https://gitlab.mbb.univ-montp2.fr/khalid/poolcnvcomp/-/tree/master

### Subset files on scaffolds

In the step 0, the goal is simply to subset strain files on scaffolds.

### Identify distant and everted read-pairs

The step 1 aime to identify distant and everted read-pairs

```{bash}
mkdir /your-path/PoolSeq_Clec/CNV/step1
mkdir /your-path/PoolSeq_Clec/CNV/step1/data
```

To go faster, we can create one script per pool (i.e. per population) using:
```{bash}
nano step0and1_LL.sh
```

```{bash}
#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem 10G
#SBATCH -t 24:00:00
#SBATCH -e /your-path/PoolSeq_Clec/CNV/step1_LL.err
#SBATCH -o /your-path/PoolSeq_Clec/CNV/step1_LL.out
#SBATCH -J Genome_Cimex_lectularius

cd /your-path/PoolSeq_Clec/CNV/

function step1d(){

    chr=$1
    pool=$2

    # Original sequence data available on Shot Read Archive, bioproject PRJNA603262
    ficBam=/your-path/PoolSeq_Clec/CNV/step1/${pool}_${chr}_mapped.bam

    /your-path/Tools/samtools-1.10/bin/samtools view --threads 8 $ficBam | awk 'function abs(x){return ((x < 0.0) ? -x : x)}{if (abs($9) < 2000 && $5>=20) print abs($9)}' > /your-path/PoolSeq_Clec/CNV/step1/data/Chr$chr/ISIZE_distribution_c${chr}p${pool}.txt # Extract the insert size distribution ; replace awk '{ if ($5 >= 20 && $9 <= 2000) print $9 }'

    sort -n -k 1 /your-path/PoolSeq_Clec/CNV/step1/data/Chr$chr/ISIZE_distribution_c${chr}p${pool}.txt  | awk '{all[NR] = $0} END{print all[int(NR*0.99 - 0.01)]}' > /your-path/PoolSeq_Clec/CNV/step1/data/Chr$chr/ISIZE99_c${chr}p${pool}.txt # Compute the 99th percentile of the insert size distribution

    insertSizeCutoff=$(cat /your-path/PoolSeq_Clec/CNV/step1/data/Chr$chr/ISIZE99_c${chr}p${pool}.txt) # Define the variable insertSizeCutoff

    /your-path/Tools/samtools-1.10/bin/samtools view --threads 8 $ficBam | python /your-path/PoolSeq_Clec/CNV/step1/step1_findDistantInserts.py $insertSizeCutoff > /your-path/PoolSeq_Clec/CNV/step1/data/Chr$chr/c${chr}p${pool}d.txt

    awk '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { print sqrt(mean2 / NR); }' /your-path/PoolSeq_Clec/CNV/step1/data/Chr$chr/ISIZE_distribution_c${chr}p${pool}.txt  > /your-path/PoolSeq_Clec/CNV/step1/data/Chr$chr/SD_ISIZE_c${chr}p${pool}.txt # calculate SD of insert size

}

# Define function for running the script to identify everted reads:

function step1e(){

    chr=$1
    pool=$2

    ficBam=/your-path/PoolSeq_Clec/CNV/step1/${pool}_${chr}_mapped.bam

    /your-path/Tools/samtools-1.10/bin/samtools flagstat --threads 8 $ficBam | sed '9q;d'| cut -d ' ' -f 1 | cat > /your-path/PoolSeq_Clec/CNV/step1/data/Chr$chr/d${pool}_c${chr}.txt # Save read depth information for later calculation of normConst ; 9q;d replaced by 7q;d ??

    /your-path/Tools/samtools-1.10/bin/samtools view --threads 8 $ficBam | python /your-path/PoolSeq_Clec/CNV/step1/step1_findEvertedInserts.py > /your-path/PoolSeq_Clec/CNV/step1/data/Chr$chr/c${chr}p${pool}e.txt # Identify everted reads

    rm -rf $ficBam
}

for chr in $(cat /your-path/PoolSeq_Clec/CNV/list_scaffolds) # iterate across all chromosomes
do
  mkdir /your-path/PoolSeq_Clec/CNV/step1/data/Chr$chr
  cd /your-path/PoolSeq_Clec/CNV/step1/ # create subset bam file here

for pool in LL # iterate across all pools
  do

     /your-path/Tools/samtools-1.10/bin/samtools view --threads 8 -b /your-path/PoolSeq_Clec/Mapped_GCF/SEP_MAP_UNMAP/${pool}_mapped_sorted.bam "${chr}" > ${pool}_${chr}_mapped.bam

     step1d $chr $pool
     step1e $chr $pool

  done
  cd ..
done
```

### Step 2 & 3: Cluster distant and everted read-pairs within pools to identify deletions and tandem duplications

I created an R script to replace steps 2 and 3 as created by North et al. 2019, with the following modifications:
- change parameter insertSizeDiffCutoff = 4*SD on .sh file
- remove computation of normcov
- remove events with frequency <2 *in both populations* (after matching between populations)

Some of these scaffolds carried events with huge number of reads, which made the soft very slow.
Those events seemed to be quite small, as here: 

NW_019392726.1,966231,966631	NW_019392726.1,966231,966631	1711.1256	NW_019392726.1,966231,966583	1224.4709	2217279	281	3443522	281	200	486.655		Quantitative	T	5.16945 	-12.38629	everted	FP	FP

The normalized read number in pop1 (i.e. LL) = 1711.1256, in pop2 (i.e. LF) = 1224.4709, size = 966631-966231 = 400)

So we decided to filter out those events based on size >=500 before running step2 and step3, since we will filter events >=1000 afterward.
We had to do it on scaffolds NW_019392706.1, NW_019392715.1, NW_019392726.1, NW_019392782.1, NW_019392787.1, NW_019392930.1, and NW_019392955.1

Here an example for NW_019392782.1:

```{r}
# Back to step 1 !
# cd /your-path/PoolSeq_Clec/CNV/step1/data/ChrNW_019392782.1
# We have to do this onto the four strains - but here we show it only on London Lab (LL)

LL <- read.table("cNW_019392782.1pLLe.txt", fill=T)
colnames(LL)=c("readid","scaffold","sl","el","strand1","qual1","seq1","seqqual1","sr","er","strand2","qual2","seq2","seqqual2")

library(dplyr)
LL <- LL %>% 
  rowwise() %>%
  mutate(left_start = min(sl,el,er,sr))
LL <- LL %>% 
  rowwise() %>%
  mutate(right_end = max(sl,el,er,sr))

LL$size <- LL$right_end - LL$left_start
dim(LL) # 1303307
LL <- LL[LL$size>=500,] 
LL <- LL[!is.na(LL$size),] 
dim(LL) # 4639

# remove extra columns
# and dowload it
LL <- as.data.frame(LL)
LL <- LL[,c(1:14)]
write.table(LL, "cNW_019392782.1pLLe_cut.txt", quote=F, sep="\t", row.names=F, col.names=F)
```

Let's create our R file with steps 2 and 3 now:

```{bash}
nano step2and3.R
```

```{r}
#!/usr/bin/env Rscript

##################################################################################################################
# parse arguments given to the script
##################################################################################################################

args = commandArgs(trailingOnly=TRUE)

# test if there is the required number of argument: if not, return an error
if (length(args)!=10) {
  stop("Ten arguments must be supplied", call.=FALSE)
}
if (length(args)==10) {
  # file1="cNW_019392689.1pLLinv.txt" # output step1 - poolA
  # file2="cNW_019392689.1pLFinv.txt" # output step1 - poolB
  # normConstA=1.31
  # normConstB=0.81
  # insertSizeDiffCutoff_popA=250 # bp difference authorized for CLUSTERED intervals
  # insertSizeDiffCutoff_popB=250 # bp difference authorized for CLUSTERED intervals
  # distanceCutoff=500 # bp difference authorized for MATCHED intervals
  # output_file_1="LL_step2.txt"
  # output_file_2="LF_step2.txt"
  # output_file_3="LL_vs_LF_step3.txt"

  file1=args[1]
  file2=args[2]
  normConstA=as.numeric(args[3])
  normConstB=as.numeric(args[4])
  insertSizeDiffCutoff_popA=as.numeric(args[5]) 
  insertSizeDiffCutoff_popB=as.numeric(args[6]) 
  distanceCutoff=as.numeric(args[7])
  output_file_1=args[8] 
  output_file_2=args[9] 
  output_file_3=args[10] 
}

##################################################################################################################
#  import data and extract usefull information
##################################################################################################################

# Load librairies
suppressPackageStartupMessages({
  library(IRanges)
  library(igraph)
  library(dplyr)
  library(stringr)
})

Pop1=read.csv2(file1, sep="\t", header=F)
Pop2=read.csv2(file2, sep="\t", header=F)

# get coordinates

Pop1_coord=Pop1[,c(2,3,4,9,10)]
names(Pop1_coord)=c("scaffold","sl","el","sr","er")  

Pop2_coord=Pop2[,c(2,3,4,9,10)]
names(Pop2_coord)=c("scaffold","sl","el","sr","er")   

scaffold <- unique(Pop1_coord$scaffold)

##################################################################################################################
# Find overlapping ranges within a population and test whether they should be clustered
##################################################################################################################

# Goal: compute number of clustered reads
# We want to find overlapping window on extreme left start and right end

# Check that they are the respective smallest/highest values

Pop1_coord <- Pop1_coord %>% 
  rowwise() %>%
  mutate(left_start = min(sl,el,er,sr))
Pop1_coord <- Pop1_coord %>% 
  rowwise() %>%
  mutate(right_end = max(sl,el,er,sr))

Pop2_coord <- Pop2_coord %>% 
  rowwise() %>%
  mutate(left_start = min(sl,el,er,sr))
Pop2_coord <- Pop2_coord %>% 
  rowwise() %>%
  mutate(right_end = max(sl,el,er,sr))

Pop1_coord$interval=paste0(Pop1_coord$left_start, "-", Pop1_coord$right_end)
Pop2_coord$interval=paste0(Pop2_coord$left_start, "-", Pop2_coord$right_end)

# Define intervals using IRanges package
intervals_pop1=IRanges(Pop1_coord$left_start, Pop1_coord$right_end)
intervals_pop2=IRanges(Pop2_coord$left_start, Pop2_coord$right_end)

# We need a match between intervals_pop1 > to cluster them

overlap1 = findOverlapPairs(query = intervals_pop1, subject = intervals_pop1)

tab_coord_x1=as.data.frame(overlap1)
tab_coord_x1$diff1=abs(tab_coord_x1$first.start-tab_coord_x1$second.start)
tab_coord_x1$diff2=abs(tab_coord_x1$first.end-tab_coord_x1$second.end)

# rule : the difference between the positions of inserts must be < insertSizeDiffCutoff
# at least in one of the two position (in the other position we allowed insertSizeDiffCutoff*2)
tab_coord_x1$match=rep(0, length(overlap1))
tab_coord_x1$match[(tab_coord_x1$diff1<insertSizeDiffCutoff_popA & tab_coord_x1$diff2<insertSizeDiffCutoff_popA) | 
                     (tab_coord_x1$diff1<insertSizeDiffCutoff_popA & tab_coord_x1$diff2<insertSizeDiffCutoff_popA)] <- 1

matches1=data.frame(as.data.frame(overlap1), tab_coord_x1$match)

matches1$first_interval=paste0("A_",matches1$first.start, "-", matches1$first.end)
matches1$second_interval=paste0("B_",matches1$second.start, "-", matches1$second.end)
g1 <- graph.edgelist( as.matrix(matches1[matches1$tab_coord_x1.match==1,][c("first_interval", "second_interval")]) )

# redefine the positions of each interval (merging of the x matched intervals ) 
# and create ouput table for clustered intervals

liste1=split(names(igraph::clusters(g1)$membership), igraph::clusters(g1)$membership)

liste_intervals1=lapply(X = liste1, FUN=function(x){
  n=length(x)
  table=strsplit(x, split = "_")
  table=do.call(rbind.data.frame,table)
  intervals=do.call(rbind.data.frame, strsplit(as.character(table[,2]), "-"))
  table=data.frame(table, intervals)
  names(table)=c("pop", "interval", "start", "end")
  table_A=table[table$pop=="A",]
  table_B=table[table$pop=="B",]
  
  table=rbind(table_A, table_B)
  table$start=as.numeric(as.character(table$start))
  table$end=as.numeric(as.character(table$end))
  
  # define new interval positions
  new_interval=reduce(IRanges(start = min(table$start), end = max(table$end))) 
  
  # append intervals from A and from B, and the number of reads supporting each
  res=data.frame(paste(scaffold, start(new_interval), end(new_interval),sep=","), 
                 paste(x=table_A$interval, sep=" ", collapse=","),
                 nreads_A=length(table_A$interval),
                 paste(x=table_B$interval, sep=" ", collapse=","),
                 nreads_B=length(table_B$interval)
  )
  names(res)=c("interval", "intervals_A", "nreads_A", "intervals_B", "nreads_B")
  return(res)
})

liste_intervals1=do.call(rbind.data.frame, liste_intervals1)
#check that columns 4:5 are just duplicates of columns 2:3
# all(liste_intervals1$nreads_A == liste_intervals1$nreads_B)

liste_intervals1 <- liste_intervals1[,c(1:3)]

# Now, for pop2:
overlap2 = findOverlapPairs(query = intervals_pop2, subject = intervals_pop2)

tab_coord_x2=as.data.frame(overlap2)
tab_coord_x2$diff1=abs(tab_coord_x2$first.start-tab_coord_x2$second.start)
tab_coord_x2$diff2=abs(tab_coord_x2$first.end-tab_coord_x2$second.end)

tab_coord_x2$match=rep(0, length(overlap2))
tab_coord_x2$match[(tab_coord_x2$diff1<insertSizeDiffCutoff_popB & tab_coord_x2$diff2<insertSizeDiffCutoff_popB) | 
                     (tab_coord_x2$diff1<insertSizeDiffCutoff_popB & tab_coord_x2$diff2<insertSizeDiffCutoff_popB)] <- 1

matches2=data.frame(as.data.frame(overlap2), tab_coord_x2$match)

matches2$first_interval=paste0("A_",matches2$first.start, "-", matches2$first.end)
matches2$second_interval=paste0("B_",matches2$second.start, "-", matches2$second.end)
g2 <- graph.edgelist( as.matrix(matches2[matches2$tab_coord_x2.match==1,][c("first_interval", "second_interval")]) )

# redefine the positions of each interval (merging of the x matched intervals) 
# and create ouput table for clustered intervals

liste2=split(names(igraph::clusters(g2)$membership), igraph::clusters(g2)$membership)

liste_intervals2=lapply(X = liste2, FUN=function(x){
  n=length(x)
  table=strsplit(x, split = "_")
  table=do.call(rbind.data.frame,table)
  intervals=do.call(rbind.data.frame, strsplit(as.character(table[,2]), "-"))
  table=data.frame(table, intervals)
  names(table)=c("pop", "interval", "start", "end")
  table_A=table[table$pop=="A",]
  table_B=table[table$pop=="B",]
  
  table=rbind(table_A, table_B)
  table$start=as.numeric(as.character(table$start))
  table$end=as.numeric(as.character(table$end))
  
  # define new interval positions
  new_interval=reduce(IRanges(start = min(table$start), end = max(table$end))) 

  # append intervals from A and from B, and the number of reads supporting each
  res=data.frame(paste(scaffold, start(new_interval), end(new_interval),sep=","), 
                 paste(x=table_A$interval, sep=" ", collapse=","),
                 nreads_A=length(table_A$interval),
                 paste(x=table_B$interval, sep=" ", collapse=","),
                 nreads_B=length(table_B$interval)
  )
  names(res)=c("interval", "intervals_A", "nreads_A", "intervals_B", "nreads_B")
  return(res)
})

liste_intervals2=do.call(rbind.data.frame, liste_intervals2)
#check that columns 4:5 are just duplicates of columns 2:3
#all(liste_intervals2$nreads_A == liste_intervals2$nreads_B)

liste_intervals2 <- liste_intervals2[,c(1:3)]

print(paste("Congrats!! Clustering (step2) is done on scaffold", scaffold))

##################################################################################################################
# Filter on the number of reads supporting the event and transform to clean tables
##################################################################################################################

threshold_reads=1 # useless here, but you can set it to higher value
intervals1_clean=liste_intervals1[liste_intervals1$nreads_A>=threshold_reads,]
intervals2_clean=liste_intervals2[liste_intervals2$nreads_A>=threshold_reads,]

if(nrow(intervals1_clean)!=0){
  intervals1_clean[c("scaffold","start","end")] <- str_split_fixed(intervals1_clean$interval, ",", 3)
  
  intervals1_clean <- intervals1_clean[,-c(1,2)]
  colnames(intervals1_clean) <- c("nreads","scaffold","start","end")
  
  intervals1_clean$start=as.numeric(as.character(intervals1_clean$start))
  intervals1_clean$end=as.numeric(as.character(intervals1_clean$end))
  
  write.table(x = intervals1_clean, file = output_file_1, quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
}else{
  print("Sadly, no cluster passed the threshold >=2 reads in 1st population.")
}

# Pop2 

if(nrow(intervals1_clean)!=0){
  intervals2_clean[c("scaffold","start","end")] <- str_split_fixed(intervals2_clean$interval, ",", 3)
  
  intervals2_clean <- intervals2_clean[,-c(1,2)]
  colnames(intervals2_clean) <- c("nreads","scaffold","start","end")
  
  intervals2_clean$start=as.numeric(as.character(intervals2_clean$start))
  intervals2_clean$end=as.numeric(as.character(intervals2_clean$end))
}else{
  print("Sadly, no cluster passed the threshold >=2 reads in 2nd population.")
}

write.table(x = intervals2_clean, file = output_file_2, quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")

##################################################################################################################
# Find overlapping ranges between populations and test whether they should be matched
##################################################################################################################

# Goal: compute number of clustered reads
# We want to find overlapping window on start and end

if(nrow(intervals1_clean)!=0){
  intervals1_clean$interval=paste0(intervals1_clean$start, "-", intervals1_clean$end)
}

if(nrow(intervals2_clean)!=0){
  intervals2_clean$interval=paste0(intervals2_clean$start, "-", intervals2_clean$end)
}

# Define intervals using IRanges package
intervals_pop1=IRanges(intervals1_clean$start, intervals1_clean$end)
intervals_pop2=IRanges(intervals2_clean$start, intervals2_clean$end)

# We need to match intervals_pop1 and intervals_pop2 

overlap = findOverlapPairs(query = intervals_pop1, subject = intervals_pop2)

tab_coord_x=as.data.frame(overlap)
tab_coord_x$diff1=abs(tab_coord_x$first.start-tab_coord_x$second.start)
tab_coord_x$diff2=abs(tab_coord_x$first.end-tab_coord_x$second.end)

tab_coord_x$match=rep(0, length(overlap))
tab_coord_x$match[(tab_coord_x$diff1<distanceCutoff & tab_coord_x$diff2<distanceCutoff)] <- 1

matches=data.frame(as.data.frame(overlap), tab_coord_x$match)

if(nrow(matches)!=0){
  matches$first_interval=paste0("A_",matches$first.start, "-", matches$first.end)
  matches$second_interval=paste0("B_",matches$second.start, "-", matches$second.end)
  g <- graph.edgelist( as.matrix(matches[matches$tab_coord_x.match==1,][c("first_interval", "second_interval")]) )
}

# redefine the positions of each interval (merging of the x matched intervals ) 
# and create ouput table for clustered intervals

if(exists("g")){ 
  if(clusters(g)$no!=0){
    
    sub_pop1 <- intervals1_clean[,c("scaffold", "start", "end", "nreads")]
    sub_pop2 <- intervals2_clean[,c("scaffold", "start", "end", "nreads")]
    
    liste=split(names(igraph::clusters(g)$membership), igraph::clusters(g)$membership)
    liste_intervals=lapply(X = liste, FUN=function(x){
      n=length(x)
      table=strsplit(x, split = "_")
      table=do.call(rbind.data.frame,table)
      intervals=do.call(rbind.data.frame, strsplit(as.character(table[,2]), "-"))
      table=data.frame(table, intervals)
      names(table)=c("pop", "interval", "start", "end")
      table_A=table[table$pop=="A",]
      table_B=table[table$pop=="B",]
    
      # add read numbers supporting each interval
      table_A=merge(x = table_A, y = sub_pop1, all.x=TRUE, all.y=FALSE)
      table_B=merge(x = table_B, y = sub_pop2, all.x=TRUE, all.y=FALSE)
      # sort on read numbers
      table_A=table_A[order(as.numeric(as.character(table_A$nreads)), decreasing = TRUE),]
      table_B=table_B[order(as.numeric(as.character(table_B$nreads)), decreasing = TRUE),]
    
      table=rbind(table_A, table_B)
      table$start=as.numeric(as.character(table$start))
      table$end=as.numeric(as.character(table$end))
    
      # sort on read numbers
      table=table[order(table$nreads, decreasing = TRUE),]
    
      # define new interval positions
      new_interval=reduce(IRanges(start = table$start, end = table$end))
    
      # append intervals from A and from B, and the number of reads supporting each
      res=data.frame(paste(table_A$scaffold[1], start(new_interval), end(new_interval),sep=","), 
                   paste(x=table_A$interval, sep=" ", collapse=","),
                   sum(table_A$nreads),
                   paste(x=table_B$interval, sep=" ", collapse=","),
                   sum(table_B$nreads)
      )
      names(res)=c("interval", "intervals_A", "nreads_A", "intervals_B", "nreads_B")
      return(res)
    })
    liste_intervals=do.call(rbind.data.frame, liste_intervals)
  
    intervals_clean <- liste_intervals
  
    intervals_clean[c("scaffold","start","end")] <- str_split_fixed(intervals_clean$interval, ",", 3)
  
    intervals_clean <- intervals_clean[,-1]
    colnames(intervals_clean) <- c("reads_A","nreads_A","reads_B","nreads_B","scaffold","start","end")
    intervals_clean$interval=paste0(intervals_clean$start, "-", intervals_clean$end)
  }
}

##################################################################################################################
# now we need to add the intervals that are specific to one population (the previous ones were shared by both)
##################################################################################################################

if(exists("intervals_clean")){
  matched_intervals_A=unlist(strsplit(as.character(intervals_clean$reads_A), split = ","))
  matched_intervals_B=unlist(strsplit(as.character(intervals_clean$reads_B), split = ","))
}

if(exists("matched_intervals_A")){
  if(!is.null(matched_intervals_A)){
    intervals1_clean_specific=intervals1_clean[-which(intervals1_clean$interval %in% matched_intervals_A),]
    }else{
      intervals1_clean_specific=intervals1_clean
    }
}

if(exists("matched_intervals_A")){
  if(!is.null(matched_intervals_B)){
    intervals2_clean_specific=intervals2_clean[-which(intervals2_clean$interval %in% matched_intervals_B),]
    }else{
      intervals2_clean_specific=intervals2_clean
    }
}

if(exists("intervals1_clean_specific")){
  if(nrow(intervals1_clean_specific)!=0){
    Pop1_coord_specific=data.frame(paste(intervals1_clean_specific$scaffold, intervals1_clean_specific$start, intervals1_clean_specific$end, sep=","),
                                 paste(intervals1_clean_specific$scaffold, intervals1_clean_specific$start, intervals1_clean_specific$end, sep=","),
                                 intervals1_clean_specific$nreads,
                                 NA,
                                 0)
    names(Pop1_coord_specific)=c("interval", "intervals_A", "nreads_A", "intervals_B", "nreads_B")
  }
}

if(exists("intervals1_clean_specific")){
  if(nrow(intervals2_clean_specific)!=0){
    Pop2_coord_specific=data.frame(paste(intervals2_clean_specific$scaffold, intervals2_clean_specific$start, intervals2_clean_specific$end, sep=","),
                                 NA,
                                 0,
                                 paste(intervals2_clean_specific$scaffold, intervals2_clean_specific$start, intervals2_clean_specific$end, sep=","),
                                 intervals2_clean_specific$nreads)
    names(Pop2_coord_specific)=c("interval", "intervals_A", "nreads_A", "intervals_B", "nreads_B")
  }
}

# combine shared and specific intervals:
if(exists("intervals_clean")){ 
  if(nrow(intervals_clean)!=0){
  All_pops_coord_specific=data.frame(paste(intervals_clean$scaffold, intervals_clean$start, intervals_clean$end, sep=","),
                                     paste(intervals_clean$scaffold, intervals_clean$start, intervals_clean$end, sep=","),
                                     intervals_clean$nreads_A,
                                     paste(intervals_clean$scaffold, intervals_clean$start, intervals_clean$end, sep=","),
                                     intervals_clean$nreads_B)
  names(All_pops_coord_specific)=c("interval", "intervals_A", "nreads_A", "intervals_B", "nreads_B")
  }
}

output <- as.data.frame(rbind(if(exists("All_pops_coord_specific")) All_pops_coord_specific,
                              if(exists("Pop1_coord_specific")) Pop1_coord_specific,
                              if(exists("Pop2_coord_specific")) Pop2_coord_specific))

if(exists("output")){ 
  if(nrow(output)!=0){
    output <- output[output$nreads_A>=2 | output$nreads_B>=2,]
  }
}

# correct number of reads to take into account differences in coverage between samples
# output$nreads_A= output$nreads_A*normConstA
# output$nreads_B= output$nreads_B*normConstB

dim(output)[1]

if(dim(output)[1]==0){ 
  print("Sadly, output is empty due to no clustered (and therefore no matched) reads")
}

print(paste("Congrats!! Matching events between populations (step3) is done on scaffold", scaffold))

##################################################################################################################
# write to disk
##################################################################################################################

write.table(x = output, file = output_file_3, quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")
```

Then, run it with the following modified command:

```{R}
step2and3.R c${chr}p${poolA}e_cut.txt c${chr}p${poolB}e_cut.txt $normConstA $normConstB $insertSizeDiffCutoff_popA $insertSizeDiffCutoff_popB $distanceCutoff c${chr}P${poolA}_everted_cluster.tsv c${chr}P${poolB}_everted_cluster.tsv c${chr}P${poolA}vsP${poolB}_everted.tsv
```

You can define normConstA, normConstB, insertSizeDiffCutoff_popA, insertSizeDiffCutoff_popB and distanceCutoff with the same parameter as we do - see in "parse arguments given to the script".

### Step 4: Calculate relative read depth differences

```{bash}
mkdir /your-path/PoolSeq_Clec/CNV/step4
mkdir /your-path/PoolSeq_Clec/CNV/step4/data
mkdir /your-path/PoolSeq_Clec/CNV/step4/data/comparison1
```

step4_comp1.sh

```{bash}
#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem 1G
#SBATCH -t 10:00:00
#SBATCH -e /your-path/PoolSeq_Clec/CNV/step4_comp1.err
#SBATCH -o /your-path/PoolSeq_Clec/CNV/step4_comp1.out
#SBATCH -J Genome_Cimex_lectularius

#Implement step4 (d and e) for comparisons 10 and 11 (the two control comparisons)

function step4() {
chr=$1
comparison=$2

# Read in the masked-regions file for the relevant chromosome (e.g repeats)
MaskedRegionsFile=/your-path/PoolSeq_Clec/CNV/bed_files/Cimex_lectularius_TE_Chr${chr}.bed

# Define comparisons
if [[ ${comparison} -eq 1 ]]; then # Focal comparison 1
  poolA=LL
  poolB=LF
elif [[ ${comparison} -eq 2 ]]; then # Control comparison 1
  poolA=LL
  poolB=GL 
elif [[ ${comparison} -eq 3 ]]; then # Focal comparison 2
  poolA=LL
  poolB=SF
elif [[ ${comparison} -eq 4 ]]; then # Control comparison 2
  poolA=LF
  poolB=SF
elif [[ ${comparison} -eq 5 ]]; then
  poolA=GL
  poolB=LF
elif [[ ${comparison} -eq 6 ]]; then
  poolA=GL
  poolB=SF
# 10 = also control comparison 1
# 11 = also control comparison 2
fi

# Read in the relevant .bam file
ficBamA=/your-path/PoolSeq_Clec/CNV/step4/${poolA}_${chr}_mapped.bam
ficBamB=/your-path/PoolSeq_Clec/CNV/step4/${poolB}_${chr}_mapped.bam

# distant read-pair clusters
python /your-path/PoolSeq_Clec/CNV/step4/step4_countReadPairsInCNV-KH-pySam.py \
/your-path/PoolSeq_Clec/CNV/step3/data/comparison${comparison}/c${chr}P${poolA}vsP${poolB}_distant_clean.tsv \
  $ficBamA $MaskedRegionsFile > /your-path/PoolSeq_Clec/CNV/step4/data/comparison${comparison}/c${chr}_Aonly_distant.tsv # this is an intermediate file

python /your-path/PoolSeq_Clec/CNV/step4/step4_countReadPairsInCNV-KH-pySam.py \
/your-path/PoolSeq_Clec/CNV/step4/data/comparison${comparison}/c${chr}_Aonly_distant.tsv \
  $ficBamB $MaskedRegionsFile  > /your-path/PoolSeq_Clec/CNV/step4/data/comparison${comparison}/c${chr}_P${poolA}vsP${poolB}d_comp${comparison}_distant.tsv

# everted read-pair clusters only
python2 /your-path/PoolSeq_Clec/CNV/step4/step4_countReadPairsInCNV-KH-pySam.py \
/your-path/PoolSeq_Clec/CNV/step2and3/data/comparison${comparison}/c${chr}P${poolA}vsP${poolB}_everted.tsv \
  $ficBamA $MaskedRegionsFile > /your-path/PoolSeq_Clec/CNV/step4/data/comparison${comparison}/c${chr}_Aonly_everted.tsv # this is an intermediate file

python2 /your-path/PoolSeq_Clec/CNV/step4/step4_countReadPairsInCNV-KH-pySam.py \
/your-path/PoolSeq_Clec/CNV/step4/data/comparison${comparison}/c${chr}_Aonly_everted.tsv \
  $ficBamB $MaskedRegionsFile  > /your-path/PoolSeq_Clec/CNV/step4/data/comparison${comparison}/c${chr}_P${poolA}vsP${poolB}d_comp${comparison}_everted.tsv

# remove the relevant data at the end of this part of the loop
rm -f   $ficBamA
rm -f   $ficBamB
}
# comparison 1 (Our focal comparison 1), for example
for chr in $(cat /your-path/PoolSeq_Clec/CNV/list_scaffolds)
do
  for comparison in 1 # only comparison 1 here
  do
    /your-path/Tools/samtools-1.10/bin/samtools view --threads 8 -b /your-path/PoolSeq_Clec/Mapped_GCF/SEP_MAP_UNMAP/LL_mapped_sorted.bam "${chr}" > /your-path/PoolSeq_Clec/CNV/step4/LL_${chr}_mapped.bam
    
    /your-path/Tools/samtools-1.10/bin/samtools index /your-path/PoolSeq_Clec/CNV/step4/LL_${chr}_mapped.bam
    
    /your-path/Tools/samtools-1.10/bin/samtools view --threads 8 -b /your-path/PoolSeq_Clec/Mapped_GCF/SEP_MAP_UNMAP/LF_mapped_sorted.bam "${chr}" > /your-path/PoolSeq_Clec/CNV/step4/LF_${chr}_mapped.bam
    
    /your-path/Tools/samtools-1.10/bin/samtools index /your-path/PoolSeq_Clec/CNV/step4/LF_${chr}_mapped.bam
    
    step4 $chr $comparison
  done
done
```


### Step 4.5: Append additional allele frequency information to each duplication and deletion event

We kept the same comparison name as in step4 > so we replace "altname" with "comparison" in for example "~/step4/data/comparison${altname}/c${chr}_P${poolA}vsP${poolB}d_comp${altname}_distant.tsv"

So first, we want to compute upper 95% threshold and lower 5% threshold for the *expected difference in
read depth*, based on empirical sampling of the data:

```{bash}
thresholdsfile="/your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/C"${comparison}"/Up90_allChr_C"${comparison}".tsv"
thresholdsfile="/your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/C"${comparison}"/Lo10_allChr_C"${comparison}".tsv"
```

We checked step2and3 output:

For each comparison, we have an everted and a distant file with columns:
V3 = The (corrected) number of read pairs supporting the event in pool 1 = LL (0 if only found in pool 2)
V5 = The (corrected) number of read pairs supporting the event in pool 2 = LF (0 if only found in pool 2)

quantile(cNW_019392721.1PPLLvsPLF_everted$V3, probs=seq(0,1,0.01))

To create random intervals per size class, we used bedtools, with the file GCF_000648675.2_Clec_2.1_genomic.txt which contains the name of the scaffolds and their size:

```{bash}
mkdir random_window
cd /your-path/PoolSeq_Clec/CNV

/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 500 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-500.txt

/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 1000 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-1kb.txt

/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 2000 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-2kb.txt

/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 5000 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-5kb.txt

/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 8000 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-8kb.txt

/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 10000 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-10kb.txt

/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 20000 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-20kb.txt

/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 30000 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-30kb.txt

/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 40000 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-40kb.txt

/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 50000 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-50kb.txt

# Added the 05/09/22
/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 100000 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-100kb.txt

/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 500000 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-500kb.txt

/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 1000000 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-1Mb.txt

/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 1500000 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-1.5Mb.txt

/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 2000000 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-2Mb.txt

/your-path/Tools/bedtools2/bin/bedtools random -seed 123456 -n 1000 -l 2500000 -g GCF_000648675.2_Clec_2.1_genomic.txt > /your-path/PoolSeq_Clec/CNV/random-wind/random-windows-2.5Mb.txt
```

Here, we created intervals of length 1500bp:

NW_019392727.1	514869	516369	1	1500	-
NW_019392833.1	711623	713123	2	1500	+

We moved files in diff. directories:

```{bash}
for SizeClass in $(cat /your-path/PoolSeq_Clec/CNV/ListeSize)
do 
mkdir "${SizeClass}"
done

for SizeClass in $(cat /your-path/PoolSeq_Clec/CNV/ListeSize);
do
for file in ./random-windows-"${SizeClass}".txt;
do 
mv "${file}" "${SizeClass}"
done
done
```

We need one file per scaffold : 

```{bash}
cd /your-path/PoolSeq_Clec/CNV/random-wind/

for SizeClass in $(cat /your-path/PoolSeq_Clec/CNV/ListeSize);
do
cd "${SizeClass}"
for file in ./random-windows-*.txt;
do 
awk -F'\t' '{print > $1".txt"}' "${file}"
cd ..
done
done
```

But we want it to be formatted like the output of step2and3:

NW_019394230.1,43996,44176	NW_019394230.1,43996,44176	2	NW_019394230.1,43996,44176	2

```{bash}
for SizeClass in $(cat /your-path/PoolSeq_Clec/CNV/ListeSize);
do
cd $SizeClass
for file in ./NW_*.txt;
do 
awk '{print $1","$2","$3"\t"$1","$2","$3"\t"$4"\t"$1","$2","$3"\t"$4}' $file > $file"2"
done
cd ..
done
```

Faire tourner step4 sur :
/your-path/PoolSeq_Clec/CNV/random-wind/${SizeClass}/${chr}.txt

```{bash}
#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem 1G
#SBATCH -t 10:00:00
#SBATCH -e /your-path/PoolSeq_Clec/CNV/step4_randomwind_comp1.err
#SBATCH -o /your-path/PoolSeq_Clec/CNV/step4_randomwind_comp1.out
#SBATCH -J Genome_Cimex_lectularius

#Implement step4 (d and e) for comparisons 10 and 11 (the two control comparisons)

function step4() {
  chr=$1
  comparison=$2
  
  # Read in the masked-regions file for the relevant chromosome (e.g repeats)
  MaskedRegionsFile=/your-path/PoolSeq_Clec/CNV/bed_files/Cimex_lectularius_TE_Chr${chr}.bed
  
  # Define comparisons
  if [[ ${comparison} -eq 1 ]]; then # Focal comparison 1
  poolA=LL
  poolB=LF
  elif [[ ${comparison} -eq 2 ]]; then # Control comparison 1
  poolA=LL
  poolB=GL 
  # 10 = also control comparison 1
  # 11 = also control comparison 2
  fi
  
  # Read in the relevant .bam file
  ficBamA=/your-path/PoolSeq_Clec/CNV/random-wind/${SizeClass}/${poolA}_${chr}_mapped.bam
  ficBamB=/your-path/PoolSeq_Clec/CNV/random-wind/${SizeClass}/${poolB}_${chr}_mapped.bam
  
  # everted read-pair clusters only
  python2 /your-path/PoolSeq_Clec/CNV/step4/step4_countReadPairsInCNV-KH-pySam.py \
  /your-path/PoolSeq_Clec/CNV/random-wind/${SizeClass}/${chr}.txt2 \
  $ficBamA $MaskedRegionsFile > /your-path/PoolSeq_Clec/CNV/random-wind/${SizeClass}/${chr}_Aonly_random.tsv # this is an intermediate file
  
  python2 /your-path/PoolSeq_Clec/CNV/step4/step4_countReadPairsInCNV-KH-pySam.py \
  /your-path/PoolSeq_Clec/CNV/random-wind/${SizeClass}/${chr}_Aonly_random.tsv \
  $ficBamB $MaskedRegionsFile  > /your-path/PoolSeq_Clec/CNV/random-wind/${SizeClass}/${chr}_${poolA}vsP${poolB}_random.tsv
  
  # remove the relevant data at the end of this part of the loop
  rm -f   $ficBamA
  rm -f   $ficBamB
  rm -f   ${ficBamA}.bai
  rm -f   ${ficBamB}.bai
}
# comparison 1 (Our focal comparison 1), for example
for SizeClass in $(cat /your-path/PoolSeq_Clec/CNV/ListeSize)
do
for chr in $(cat /your-path/PoolSeq_Clec/CNV/list_scaffolds)
do
for comparison in 1 # only comparison 1 here
do
/your-path/Tools/samtools-1.10/bin/samtools view --threads 8 -b /your-path/PoolSeq_Clec/Mapped_GCF/SEP_MAP_UNMAP/LL_mapped_sorted.bam "${chr}" > /your-path/PoolSeq_Clec/CNV/random-wind/${SizeClass}/LL_${chr}_mapped.bam

/your-path/Tools/samtools-1.10/bin/samtools index /your-path/PoolSeq_Clec/CNV/random-wind/${SizeClass}/LL_${chr}_mapped.bam

/your-path/Tools/samtools-1.10/bin/samtools view --threads 8 -b /your-path/PoolSeq_Clec/Mapped_GCF/SEP_MAP_UNMAP/LF_mapped_sorted.bam "${chr}" > /your-path/PoolSeq_Clec/CNV/random-wind/${SizeClass}/LF_${chr}_mapped.bam

/your-path/Tools/samtools-1.10/bin/samtools index /your-path/PoolSeq_Clec/CNV/random-wind/${SizeClass}/LF_${chr}_mapped.bam

step4 $chr $comparison
done
done
done
```

Then, merge all files per scaffold into one:

```{bash}
cat NW_019*.1_PLLvsPLF_random.tsv > PLLvsPLF_random.txt
```

We then compute the read depth ratio as in step4.5 (12th field):
(We add conditions in (1) and (2) to avoid the case "divided by 0" which results in fatal error)

(1) calculate the per-nucleotide read depth in pool A if the CNV was present in A:

    Wsize=$7; depth=$6; if ($7 != 0) print (depth/Wsize)*100 ; else if ($7 == 0) print "0"
    (depth/Wsize)*100 > $1
    
(2) calculate the per-nucleotide read depth in pool B if the CNV was present in B

    Wsize=$7; depth=$8;  if ($7 != 0) print (depth/Wsize)*100 ; else if ($7 == 0) print "0"
    (depth/Wsize)*100 > $2
  
(3) compute the ratio A/B

    if ($2 != 0) print $1/$2; else if ($2 == 0) print "inf"  

Do this computation on R:

```{r}
test <- read.table("PLLvsPLF_random.txt")
test$ratio <- (test$V6/test$V7)/(test$V8/test$V9)
quantile(test$ratio,c(0.25,0.75), na.rm=T)
q()
n
```

```{bash}
cd /your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95
```

```{r}
LL1kb <- read.table("out-cov-LL-1kb.txt")
LF1kb <- read.table("out-cov-LF-1kb.txt")

colnames(LL1kb) <- c("rname","startpos","endpos","numreads","covbases")
colnames(LF1kb) <- c("rname","startpos","endpos","numreads","covbases")

LL1kb$ratio <- (LL1kb$numreads/LL1kb$covbases*100)/(LF1kb$numreads/LF1kb$covbases*100)
quantile(LL1kb$ratio, na.rm=T, prob=seq(0,1,0.05))

write.table(LL1kb, "LL1kb.txt", quote=F)
```

  V1             V2     V3     V4      ratio
1 NW_019392727.1 514869 515869 13.5235 -2.3786
2 NW_019392833.1 711623 712623 14.4935 -5.6963
3 NW_019392920.1 148931 149931 15.7602 -1.1189
4 NW_019392813.1 181956 182956 14.7433 -7.989
5 NW_019392737.1 1108705 1109705 16.7043 -12.1468
6 NW_019392689.1 912537 913537 12.6234 -3.4615

Or we can define the 5% and 95% quantiles with:

```{bash}
sort -n -k 7 LL30kb.txt  | awk '{all[NR] = $7} END{print all[int(NR*0.95)]}' > q95_30kb
sort -n -k 7 LL30kb.txt  | awk '{all[NR] = $7} END{print all[int(NR*0.05)]}' > q5_30kb
```

We can also define the 10% and 90% quantiles :

```{bash}
sort -n -k 7 LL30kb.txt  | awk '{all[NR] = $7} END{print all[int(NR*0.90)]}' > q90_30kb
sort -n -k 7 LL30kb.txt  | awk '{all[NR] = $7} END{print all[int(NR*0.10)]}' > q10_30kb
```

To do this over all files:

```{bash}
for SizeClass in $(cat /your-path/PoolSeq_Clec/CNV/ListeSize)
do
/your-path/Tools/R-4.0.5/bin/Rscript /your-path/PoolSeq_Clec/CNV/compute-quantiles-random-windows.R /your-path/PoolSeq_Clec/CNV/result-cov-random/random-windows-"$SizeClass".txt_covLL_ok.txt /your-path/PoolSeq_Clec/CNV/result-cov-random/random-windows-"$SizeClass".txt_covLF_ok.txt /your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/quantiles_"$SizeClass".txt
done
```

Then, we can manually create a file with all categories:

nano C1_PERCENTILES.txt

size    q5        q10       q90       q95
500     0.5000000 0.5789474 1.2178744 1.3539848 
1000    0.5691101 0.6234683 1.0826484 1.1762483
2000    0.6414445 0.6779114 0.9764501 1.0391094 
5000    0.7016958 0.7308842 0.9252077 0.9553209 
8000    0.7232878 0.7442450 0.9002137 0.9256592 
10000   0.7322942 0.7519393 0.8899569 0.9110690 
20000   0.7583027 0.7717086 0.8676708 0.8816779 
30000   0.7669041 0.7790163 0.8595701 0.8713869 
40000   0.7704131 0.7829543 0.8547581 0.8661766 
50000   0.7780854 0.7889317 0.8509359 0.8622906 
100000    0.7873172 0.7943371 0.8390067 0.8474167 
500000    0.8007618 0.8054869 0.8269708 0.8303853 
1000000   0.8056730 0.8079778 0.8246806 0.8272439 
1500000   0.8057133 0.8080424 0.8223257 0.8245331 
2000000   0.8068726 0.8093561 0.8209310 0.8225129 
2500000   0.8053877 0.8087843 0.8197390 0.8212512 

Create threshold files:

```{bash}
awk '{print $1,$2}' C1_PERCENTILES.txt > /your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/C1/Lo5_allChr_C1.tsv

awk '{print $1,$3}' C1_PERCENTILES.txt > /your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/C1/Lo10_allChr_C1.tsv

awk '{print $1,$4}' C1_PERCENTILES.txt > /your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/C1/Up90_allChr_C1.tsv

awk '{print $1,$5}' C1_PERCENTILES.txt > /your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/C1/Up95_allChr_C1.tsv
```

> We created several directory

```{bash}
mkdir /your-path/PoolSeq_Clec/CNV/step4pt5
mkdir /your-path/PoolSeq_Clec/CNV/step4pt5/data
mkdir /your-path/PoolSeq_Clec/CNV/step4pt5/data/comparison1
mkdir /your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/
mkdir /your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/C1
```

> We also added the line in the code about MNSI which was missing.

```{bash}
#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem 1G
#SBATCH -t 10:00:00
#SBATCH -e /your-path/PoolSeq_Clec/CNV/step4pt5_comp1.err
#SBATCH -o /your-path/PoolSeq_Clec/CNV/step4pt5_comp1.out
#SBATCH -J Genome_Cimex_lectularius

function step4pt5(){

    chr=$1
    comparison=$2

    # Define the pools relevant to the specified comparison

    if [[ ${comparison} -eq 1 ]]; then # Focal comparison 1
        poolA=LL
        poolB=LF
    elif [[ ${comparison} -eq 2 ]]; then # Focal comparison 2
        poolA=LL
        poolB=GL
    elif [[ ${comparison} -eq 3 ]]; then # Control comparison 1
        poolA=LL
        poolB=SF
    elif [[ ${comparison} -eq 4 ]]; then # Control comparison 2
        poolA=LF
        poolB=SF
    elif [[ ${comparison} -eq 5 ]]; then
        poolA=GL
        poolB=LF
    elif [[ ${comparison} -eq 6 ]]; then
        poolA=GL
        poolB=SF
    fi
    
    ### EVERTED READS (DUPLICATIONS)

    ## Assign size window categories to each CNV

    cut -f7 /your-path/PoolSeq_Clec/CNV/step4/data/comparison${comparison}/c${chr}_P${poolA}vsP${poolB}d_comp${comparison}_everted.tsv | awk -F $'\t' 'BEGIN {OFS = FS} \
     {if ($1 < 750) print "500"; \
     else if ($1 >= 750 && $1 < 1500) print "1000"; \
     else if ($1 >= 1500 && $1 < 3500) print "2000"; \
     else if ($1 >= 3500 && $1 < 6500) print "5000"; \
     else if ($1 >= 6500 && $1 < 9000) print "8000";\
     else if ($1 >= 9000 && $1 < 15000) print "10000";\
     else if ($1 >= 15000 && $1 < 25000) print "20000";\
     else if ($1 >= 25000 && $1 < 35000) print "30000";\
     else if ($1 >= 35000 && $1 < 45000) print "40000";\
     else if ($1 >= 45000 && $1 < 75000) print "50000";\
     else if ($1 >= 75000 && $1 < 300000) print "100000";\
     else if ($1 >= 300000 && $1 < 750000) print "500000";\
     else if ($1 >= 750000 && $1 < 1250000) print "1000000";\
     else if ($1 >= 1250000 && $1 < 1750000) print "1500000";\
     else if ($1 >= 1750000 && $1 < 2250000) print "2000000";\
     else if ($1 >= 2250000) print "2500000"}' > windowSizes_Chr${chr}_Pools${poolA}vs${poolB}_e.tsv

    ## Append the 10th field with window category information
    paste /your-path/PoolSeq_Clec/CNV/step4/data/comparison${comparison}/c${chr}_P${poolA}vsP${poolB}d_comp${comparison}_everted.tsv \
    windowSizes_Chr${chr}_Pools${poolA}vs${poolB}_e.tsv > S4_Window_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv

    ## Append the 11th field with the difference in number of supporting inserts (DiffSuppIns):

    awk -F $'\t' 'BEGIN {OFS = FS} {print $3-$5 }' S4_Window_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv > DSI_e.tsv
    paste S4_Window_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv DSI_e.tsv > S4_DSI_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv

     ## Calculate the read depth ratio and append the 12th field with this value
     
    # calculate the per-nucleotide read depth in pool A if the CNV was present in A
    awk -F $'\t' 'BEGIN {OFS = FS}{Wsize=$7; depth=$6; if ($7 != 0) print (depth/Wsize)*100 ; else if ($7 == 0) print "0"}' S4_DSI_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv  > PerNtDepthA_e.tsv
    # calculate the per-nucleotide read depth in pool B if the CNV was present in B
    awk -F $'\t' 'BEGIN {OFS = FS}{Wsize=$7; depth=$8; if ($7 != 0) print (depth/Wsize)*100 ; else if ($7 == 0) print "0"}' S4_DSI_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv  > PerNtDepthB_e.tsv
    paste PerNtDepthA_e.tsv PerNtDepthB_e.tsv > ObsDepth_e.tsv
    awk -F $'\t' 'BEGIN {OFS = FS}{if ($2 != 0) print $1/$2; else if ($2 == 0) print "inf" }' ObsDepth_e.tsv > ObsRDR_e.tsv
    paste S4_DSI_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv ObsRDR_e.tsv > S4_ObsRDR_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv

    ## Assign a status (qualitative difference or quantitative difference) to each CNV and apppend the 13th field with this status
    awk -F $'\t' 'BEGIN {OFS = FS}{if ($2 != "NA" && $4 != "NA") print "Quantitative"; else if ($2 = "NA" || $4 = "NA") print "Qualitative"}' \
    S4_ObsRDR_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv > QQ_e.tsv
    paste S4_ObsRDR_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv QQ_e.tsv > S4_QQ_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    # Qualitative/Quantitative status is now the 13th field

# For clusters of everted reads, identify those which have at least 3 supporting reads  (MNSI>3) and append the 14th field with this status
    awk -F $'\t' 'BEGIN {OFS = FS} {if ($2 != "NA" && $4 != "NA" && ($3 <= 3 || $5 <= 3)) print "F"; \
    else if ($2 != "NA" && $4 != "NA" && ($3 >= 3 || $5 >= 3)) print "T"; \
    else if ($2 != "NA" && $3 <= 3) print "F";\
    else if ($2 != "NA" && $3 >= 3) print "T"; \
    else if ($4 != "NA" && $5 <= 3) print "F";\
    else if ($4 != "NA" && $5 >= 3) print "T"}' S4_QQ_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv > MNSI_everted.tsv
    paste S4_QQ_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv MNSI_everted.tsv > S4_MNSI_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv

    ## Read in the upper 95% read depth ratio threshold and assign this value to the 15th field

    infile="S4_MNSI_everted_Chr"${chr}"_Pools"${poolA}"vs"${poolB}".tsv"
    outfile="S4_MNSI_everted_Chr"${chr}"_Pools"${poolA}"vs"${poolB}"_Upper.tsv"
    thresholdsfile="/your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/C"${comparison}"/Up90_allChr_C"${comparison}".tsv"

    awk -v outf=$outfile -F $'\t' 'BEGIN {OFS = FS}  { if (NR == FNR ) {d[$1]=$2}
    else {
    print $0, d[$10] > outf
    }}
    ' $thresholdsfile $infile

    # 16th field is lower threshold

    infile="S4_MNSI_everted_Chr"${chr}"_Pools"${poolA}"vs"${poolB}"_Upper.tsv"
    outfile="S4_MNSI_everted_Chr"${chr}"_Pools"${poolA}"vs"${poolB}"_Lower.tsv"
    thresholdsfile="/your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/C"${comparison}"/Lo10_allChr_C"${comparison}".tsv"

    awk -v outf=$outfile -F $'\t' 'BEGIN {OFS = FS} { if (NR == FNR ) {d[$1]=$2}
    else {
    print $0, d[$10] > outf
    }}
    ' $thresholdsfile $infile

    # 17th field labels Everted reads

    awk -F $'\t' 'BEGIN {OFS = FS} {if ($14 == "T" || $14 == "F") print "everted"}' S4_MNSI_everted_Chr${chr}_Pools${poolA}vs${poolB}_Lower.tsv > EVERTED.tsv
    paste S4_MNSI_everted_Chr${chr}_Pools${poolA}vs${poolB}_Lower.tsv EVERTED.tsv > S4_labeled_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv


    # 18th field labels legitimate tandem duplications that belong to a specific pool (DUP, for everted reads), legitimate deletions (DEL, for distant reads) and false positives (FP)

    # 18th field labels legitimate tandem duplications (DUP) in to a specific pool (A or B) and false positives (FP).
    # FPs are generated when the DiffSuppIns value is not in the same 'direction' as the read depth ratio.
    # If there is a duplication in pool A, there should be more everted read
    # pairs in A compared to B and a higher read depth in A compared to B (and vice versa for b).

    awk -F $'\t' 'BEGIN {OFS = FS} {if ($11 > 0 && $12 > 1) print "DUP_A"; \
    else if ($11 < 0 && $12 < 1) print "DUP_B"; \
    else print "FP"}' S4_labeled_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv > S4_CNVAB_e.txt

    paste S4_labeled_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv S4_CNVAB_e.txt > S4_CNVAB_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv

    # 19th field labels legitimate tandem duplications in general

    awk -F $'\t' 'BEGIN {OFS = FS} {if ($18 == "DUP_A" || $18 == "DUP_B") print "DUP"; else print "FP"}' S4_CNVAB_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv > DUP_FP.tsv
    paste S4_CNVAB_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv DUP_FP.tsv > S4_appended_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv

    # remove intermediate files

    rm -rf windowSizes_Chr${chr}_Pools${poolA}vs${poolB}_e.tsv
    rm -rf S4_Window_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    rm -rf DSI_e.tsv
    rm -rf S4_DSI_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    rm -rf PerNtDepthA_e.tsv
    rm -rf PerNtDepthB_e.tsv
    rm -rf ObsDepth_e.tsv
    rm -rf ObsRDR_e.tsv
    rm -rf EVERTED.tsv
    rm -rf S4_ObsRDR_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    rm -rf UpThreshold99.tsv
    rm -rf UpThreshold99.tsv
    rm -rf MNSI_everted.tsv
    rm -rf QQ_e.tsv
    rm -rf S4_QQ_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    rm -rf S4_MNSI_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    rm -rf S4_MNSI_everted_Chr${chr}_Pools${poolA}vs${poolB}_Upper.tsv
    rm -rf S4_MNSI_everted_Chr${chr}_Pools${poolA}vs${poolB}_Lower.tsv
    rm -rf S4_labeled_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    rm -rf S4_CNVAB_e.txt
    rm -rf S4_CNVAB_everted_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    rm -rf DUP_FP.tsv

    ### DISTANT READS (DELETIONS)

    ## Assign size window categories to each CNV
    cut -f7 /your-path/PoolSeq_Clec/CNV/step4/data/comparison${comparison}/c${chr}_P${poolA}vsP${poolB}d_comp${comparison}_distant.tsv | awk -F $'\t' 'BEGIN {OFS = FS} \
     {if ($1 < 750) print "500"; \
     else if ($1 >= 750 && $1 < 1500) print "1000"; \
     else if ($1 >= 1500 && $1 < 3500) print "2000"; \
     else if ($1 >= 3500 && $1 < 6500) print "5000"; \
     else if ($1 >= 6500 && $1 < 9000) print "8000";\
     else if ($1 >= 9000 && $1 < 15000) print "10000";\
     else if ($1 >= 15000 && $1 < 25000) print "20000";\
     else if ($1 >= 25000 && $1 < 35000) print "30000";\
     else if ($1 >= 35000 && $1 < 45000) print "40000";\
     else if ($1 >= 45000 && $1 < 75000) print "50000";\
     else if ($1 >= 75000 && $1 < 300000) print "100000";\
     else if ($1 >= 300000 && $1 < 750000) print "500000";\
     else if ($1 >= 750000 && $1 < 1250000) print "1000000";\
     else if ($1 >= 1250000 && $1 < 1750000) print "1500000";\
     else if ($1 >= 1750000 && $1 < 2250000) print "2000000";\
     else if ($1 >= 2250000) print "2500000"}' > windowSizes_Chr${chr}_Pools${poolA}vs${poolB}.tsv

    ## Append the 10th field with window category information
    paste /your-path/PoolSeq_Clec/CNV/step4/data/comparison${comparison}/c${chr}_P${poolA}vsP${poolB}d_comp${comparison}_distant.tsv \
    windowSizes_Chr${chr}_Pools${poolA}vs${poolB}.tsv > S4_Window_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv

     ## Append the 11th field with the difference in number of supporting inserts (DiffSuppIns):
    awk -F $'\t' 'BEGIN {OFS = FS}{print $3-$5 }' S4_Window_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv > DSI.tsv
    paste S4_Window_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv DSI.tsv > S4_DSI_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv

    ## Calculate the read depth ratio and append the 12th field with this value

    # calculate the per-nucleotide read depth in pool A if the CNV was present in A
    awk -F $'\t' 'BEGIN {OFS = FS}{Wsize=$7; depth=$6; if ($7 != 0) print (depth/Wsize)*100 ; else if ($7 == 0) print "0"}' S4_DSI_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv  > PerNtDepthA.tsv
    # calculate the per-nucleotide read depth in pool B if the CNV was present in B
    awk -F $'\t' 'BEGIN {OFS = FS}{Wsize=$7; depth=$8; if ($7 != 0) print (depth/Wsize)*100 ; else if ($7 == 0) print "0"}' S4_DSI_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv  > PerNtDepthB.tsv
    paste PerNtDepthA.tsv PerNtDepthB.tsv > ObsDepth.tsv
    awk -F $'\t' 'BEGIN {OFS = FS}{if ($2 != 0) print $1/$2; else if ($2 == 0) print "inf" }' ObsDepth.tsv > ObsRDR.tsv
    paste S4_DSI_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv ObsRDR.tsv > S4_ObsRDR_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv

    # Assign a status (qualitative difference or quantitative difference) to each CNV and apppend the 13th field with this status
    awk -F $'\t' 'BEGIN {OFS = FS} {if ($2 != "NA" && $4 != "NA") print "Quantitative"; else if ($2 = "NA" || $4 = "NA") print "Qualitative"}' \
    S4_ObsRDR_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv > QQ.tsv
    paste S4_ObsRDR_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv QQ.tsv > S4_QQ_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv


    # For clusters of distant reads, identify those which have at least 3 supporting reads  (MNSI>3) and append the 14th field with this status
    awk -F $'\t' 'BEGIN {OFS = FS} {if ($2 != "NA" && $4 != "NA" && ($3 <= 3 || $5 <= 3)) print "F"; \
    else if ($2 != "NA" && $4 != "NA" && ($3 >= 3 || $5 >= 3)) print "T"; \
    else if ($2 != "NA" && $3 <= 3) print "F";\
    else if ($2 != "NA" && $3 >= 3) print "T"; \
    else if ($4 != "NA" && $5 <= 3) print "F";\
    else if ($4 != "NA" && $5 >= 3) print "T"}' S4_QQ_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv > MNSI_distant.tsv
    paste S4_QQ_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv MNSI_distant.tsv > S4_MNSI_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv

    # Read in the upper threshold  and assign this value to the 15th field

    infile="S4_MNSI_distant_Chr"${chr}"_Pools"${poolA}"vs"${poolB}".tsv"
    outfile="S4_MNSI_distant_Chr"${chr}"_Pools"${poolA}"vs"${poolB}"_Upper.tsv"
    thresholdsfile="/your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/C"${comparison}"/Up90_allChr_C"${comparison}".tsv"

    awk -v outf=$outfile -F $'\t' 'BEGIN {OFS = FS}  { if (NR == FNR ) {d[$1]=$2}
    else {
    print $0, d[$10] > outf
    }}
    ' $thresholdsfile $infile

    # 16th field is lower threshold

    infile="S4_MNSI_distant_Chr"${chr}"_Pools"${poolA}"vs"${poolB}"_Upper.tsv"
    outfile="S4_MNSI_distant_Chr"${chr}"_Pools"${poolA}"vs"${poolB}"_Lower.tsv"
    thresholdsfile="/your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/C"${comparison}"/Lo10_allChr_C"${comparison}".tsv"

    awk -v outf=$outfile -F $'\t' 'BEGIN {OFS = FS}  { if (NR == FNR ) {d[$1]=$2}
    else {
    print $0, d[$10] > outf
    }}
    ' $thresholdsfile $infile

    # 17th field labels Distant reads

    awk -F $'\t' 'BEGIN {OFS = FS} {if ($14 == "T" || $14 == "F") print "distant"}' S4_MNSI_distant_Chr${chr}_Pools${poolA}vs${poolB}_Lower.tsv > distant.tsv
    paste S4_MNSI_distant_Chr${chr}_Pools${poolA}vs${poolB}_Lower.tsv distant.tsv > S4_labeled_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv


    # 18th field labels legitimate tandem deletions (DEL) in to a specific pool (A or B) and false positives (FP).
    # FPs are generated when the DiffSuppIns value is not in the same 'direction' as the read depth ratio.
    # If there is a deletion in pool A, there should be more distant read
    # pairs in A compared to B and a lower read depth in A compared to B (and vice versa for b).

    awk -F $'\t' 'BEGIN {OFS = FS} {if ($11 > 0 && $12 < 1) print "DEL_A"; \
    else if ($11 < 0 && $12 > 1) print "DEL_B"; \
    else print "FP"}' S4_labeled_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv > CNVAB2.tsv

    paste S4_labeled_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv CNVAB2.tsv > S4_CNVAB_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv

    # 19th field labels each as DEL (legitimate deletion) vs FP (false positive)

    awk -F $'\t' 'BEGIN {OFS = FS} {if ($18 == "DEL_A" || $18 == "DEL_B") print "DEL"; else print "FP"}' S4_CNVAB_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv > DEL_FP.tsv
    paste S4_CNVAB_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv DEL_FP.tsv > S4_appended_distant_Chr${chr}_comparison${comparison}.tsv

        # All intermediate files are removed
    rm -rf windowSizes_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    rm -rf S4_Window_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    rm -rf DSI.tsv
    rm -rf S4_DSI_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    rm -rf PerNtDepthA.tsv
    rm -rf PerNtDepthB.tsv
    rm -rf ObsDepth.tsv
    rm -rf ObsRDR.tsv
    rm -rf S4_ObsRDR_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    rm -rf UpThreshold99.tsv
    rm -rf MNSI_distant.tsv
    rm -rf distant.tsv
    rm -rf QQ.tsv
    rm -rf S4_QQ_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    rm -rf S4_MNSI_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    rm -rf S4_MNSI_distant_Chr${chr}_Pools${poolA}vs${poolB}_Upper.tsv
    rm -rf S4_MNSI_distant_Chr${chr}_Pools${poolA}vs${poolB}_Lower.tsv
    rm -rf S4_labeled_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    rm -rf S4_CNVAB_distant_Chr${chr}_Pools${poolA}vs${poolB}.tsv
    rm -rf DEL_FP.tsv
    rm -rf CNVAB2.tsv

}

## focal comparisons
# comparison 1
cd /your-path/PoolSeq_Clec/CNV/step4pt5/data/comparison1
for chr in $(cat /your-path/PoolSeq_Clec/CNV/list_scaffolds)
do
 step4pt5 $chr 1
done
cat S4_appended_distant_Chr*.tsv > S4_distant_c1.tsv # collate all DEL output
cat S4_appended_everted_Chr*.tsv > S4_everted_c1.tsv # collate all DUP output
```

We can do some analyzes on R.

Deletions :

```{r}
del <- read.table(file="S4_distant_c1.tsv", fill=T) # 65 554
colnames(del) <- c("Coords","CoordsPool1","NormReadsPool1","CoordsPool2","NormReadsPool2","ReadDepthPool1","SizePool1","ReadDepthPool2","SizePool2","SizeClass","DiffSuppIns","ReadDepthRatio","Status","MNSI","UpperThreshold","LowerThreshold","Type","TypeDel","extra")

# To correct manually...
del_not_ok <- del[del$LowerThreshold=="distant",] # 51 865

del_not_ok$extra <- del_not_ok$TypeDel
del_not_ok$TypeDel <- del_not_ok$Type
del_not_ok$Type <- del_not_ok$LowerThreshold
del_not_ok$LowerThreshold <- del_not_ok$UpperThreshold
del_not_ok$UpperThreshold <- del_not_ok$MNSI
del_not_ok$MNSI <- del_not_ok$Status
del_not_ok$Status <- del_not_ok$ReadDepthRatio

del_not_ok$ReadDepthRatio <- (del_not_ok$ReadDepthPool1/del_not_ok$SizePool1*100)/
  (del_not_ok$ReadDepthPool2/del_not_ok$SizePool2*100)

del <- del[del$LowerThreshold!="distant",] # 13 209

del <- rbind(del,del_not_ok) # 65 554 again
# All good !
rm(del_not_ok)

# We filter out mitochondrial genome
# Split dups' first column
library(stringr)
del[c('scaffold', 'start_event', 'end_event')] <- str_split_fixed(del$Coords, ',', 3)
del <- del[del$scaffold != "NC_030043.1",] #  65 512

del <- del[del$SizeClass >= 1000,] # 9 319

summary(del$TypeDel)
```

Duplications :

```{r}
dup <- read.table(file="S4_everted_c1.tsv", fill=T) # 26 808 for S4_everted_c1.tsv 
colnames(dup) <- c("Coords","CoordsPool1","NormReadsPool1","CoordsPool2","NormReadsPool2","ReadDepthPool1","SizePool1","ReadDepthPool2","SizePool2","SizeClass","DiffSuppIns","ReadDepthRatio","Status","MNSI","UpperThreshold","LowerThreshold","Type","TypeDup","extra")

library(stringr)
dup[c('scaffold', 'start_event', 'end_event')] <- str_split_fixed(dup$Coords, ',', 3)
dup <- dup[dup$scaffold != "NC_030043.1",] #  26 790

# We then filter out small events (<= 500) and re-do a barplot

dup <- dup[dup$SizeClass >= 5000,] # 3 725

summary(as.factor(dup$TypeDup)) # only dupB if we do DiffSupIns threshold
summary(dup_inv$TypeDup)
summary(as.factor(dup_save_bis$TypeDup)) # only dupB if we do DiffSupIns threshold


CNV = c("DUP_TAND","DUP_INV") # MAJ DUP_INV : 14/09
in_LL = c(13,13)
in_LF = c(2058,2871)
FP = c(1654,2238)

data<-cbind(in_LF,in_LL,FP)

barplot(t(data),beside=F,col=c("chartreuse3","darkgoldenrod1","cyan3"),ylab="Number of events",
        names=CNV,las=3,horiz=F,
        ylim=c(0,6000),xlim=c(0,6),space=0.1)
legend(x="topright", legend=c("Events in London Field","Events in London Lab","False positive"),
       cex=0.9,fill=c("chartreuse3","darkgoldenrod1","cyan3"),bty="n")     

dup_save <- merge(dup, dup_inv, by=colnames(dup), all.x=T, all.y=T) # 8 847

# We want to recover the "real" Reads number
# normConstLF <- read.table(file="normConstLF_comp1.txt")
# normConstLL <- read.table(file="normConstLL_comp1.txt")
# dup_save <- merge(dup_save, normConstLF, all.x=T, by.x="scaffold", by.y="V2")
# dup_save <- merge(dup_save, normConstLL, all.x=T, by.x="scaffold", by.y="V2")
# colnames(dup_save)[23] <- "norm_LF"
# colnames(dup_save)[24] <- "norm_LL"
# dup_save$ReadsPool1 <- dup_save$NormReadsPool1/dup_save$norm_LL
# dup_save$ReadsPool2 <- dup_save$NormReadsPool2/dup_save$norm_LF

######## SAVE HERE ######## 
library(xlsx)
write.xlsx(dup_save, file="dup_save.xls", 
           row.names = F)
########################## 
```

Then, we decided to compute a frequency for events.
We computed coverage onto the positions of event borders.

First, we extracted position of resist_final table and formated it for samtools
```{r}
library(xlsx)
dup_save <- read.xlsx(file="dup_save.xls", sheetName = "Sheet1", header = T)

dup_save$start_event <- as.numeric(as.character(dup_save$start_event))
dup_save$end_event <- as.numeric(as.character(dup_save$end_event))

dup_save$start_event_bis <- dup_save$start_event+150
dup_save$end_event_bis <- dup_save$end_event-150

dup_save$samtools_start <- paste0(dup_save$scaffold, ":", dup_save$start_event, "-",
                                  dup_save$start_event_bis)
dup_save$samtools_end <- paste0(dup_save$scaffold, ":", dup_save$end_event_bis, "-",
                                  dup_save$end_event)
dup_samtools_start <- dup_save$samtools_start
dup_samtools_end <- dup_save$samtools_end

write.table(dup_samtools_start, file="dup_samtools_start.txt", 
           row.names = F, quote=F, col.names = F)
write.table(dup_samtools_end, file="dup_samtools_end.txt", 
           row.names = F, quote=F, col.names = F)
```

Then, we simply used samtools:
We want to keep rname,startpos,endpos,numreads,covbases > -f1,2,3,4,5  

```{bash}
cd /your-path/PoolSeq_Clec/Mapped_GCF/SEP_MAP_UNMAP

while read line;
do 
/your-path/Tools/miniconda3/bin/samtools coverage -r $line LL_mapped_sorted.bam | grep ^# -v | cut -f1,2,3,4,5
done < /your-path/PoolSeq_Clec/CNV/dup_samtools_start.txt > /your-path/PoolSeq_Clec/CNV/dup_coverage_start_LL.txt

while read line;
do 
/your-path/Tools/miniconda3/bin/samtools coverage -r $line LF_mapped_sorted.bam | grep ^# -v | cut -f1,2,3,4,5
done < /your-path/PoolSeq_Clec/CNV/dup_samtools_start.txt > /your-path/PoolSeq_Clec/CNV/dup_coverage_start_LF.txt

while read line;
do 
/your-path/Tools/miniconda3/bin/samtools coverage -r $line LL_mapped_sorted.bam | grep ^# -v | cut -f1,2,3,4,5
done < /your-path/PoolSeq_Clec/CNV/dup_samtools_end.txt > /your-path/PoolSeq_Clec/CNV/dup_coverage_end_LL.txt

while read line;
do 
/your-path/Tools/miniconda3/bin/samtools coverage -r $line LF_mapped_sorted.bam | grep ^# -v | cut -f1,2,3,4,5
done < /your-path/PoolSeq_Clec/CNV/dup_samtools_end.txt > /your-path/PoolSeq_Clec/CNV/dup_coverage_end_LF.txt
```

```{r}
dup_coverage_start_LL <- read.table(file="dup_coverage_start_LL.txt")
dup_coverage_end_LL <- read.table(file="dup_coverage_end_LL.txt")
dup_coverage_start_LF <- read.table(file="dup_coverage_start_LF.txt")
dup_coverage_end_LF <- read.table(file="dup_coverage_end_LF.txt")

colnames(dup_coverage_start_LL) <- c("rname","startpos","endpos","numreads_start_LL","covbases_start_LL")
colnames(dup_coverage_end_LL) <- c("rname","startpos","endpos","numreads_end_LL","covbases_end_LL")
colnames(dup_coverage_start_LF) <- c("rname","startpos","endpos","numreads_start_LF","covbases_start_LF")
colnames(dup_coverage_end_LF) <- c("rname","startpos","endpos","numreads_end_LF","covbases_end_LF")

dup_save_bis <- merge(dup_save, dup_coverage_start_LL, by.x=c("scaffold","start_event","start_event_bis"),
                      by.y=c("rname","startpos","endpos"), all.x=T, all.y=F)
dup_save_bis <- merge(dup_save_bis, dup_coverage_end_LL, by.x=c("scaffold","end_event_bis","end_event"),
                      by.y=c("rname","startpos","endpos"), all.x=T, all.y=F)
dup_save_bis <- merge(dup_save_bis, dup_coverage_start_LF,
                      by.x=c("scaffold","start_event","start_event_bis"),
                      by.y=c("rname","startpos","endpos"), all.x=T, all.y=F)
dup_save_bis <- merge(dup_save_bis, dup_coverage_end_LF, by.x=c("scaffold","end_event_bis","end_event"),
                      by.y=c("rname","startpos","endpos"), all.x=T, all.y=F)
library(dplyr)
dup_save_bis <- distinct(dup_save_bis)

dup_save_bis$cov_junction_LL <- (dup_save_bis$numreads_start_LL+dup_save_bis$numreads_end_LL)/2
dup_save_bis$cov_junction_LF <- (dup_save_bis$numreads_start_LF+dup_save_bis$numreads_end_LF)/2

dup_save_bis$freqLL <- dup_save_bis$NormReadsPool1/dup_save_bis$cov_junction_LL # change with non-norm
dup_save_bis$freqLF <- dup_save_bis$NormReadsPool2/dup_save_bis$cov_junction_LF # change with non-norm

summary(dup_save_bis$freqLL)
summary(dup_save_bis$freqLF)

# SECOND METHOD - NOT USING SAMTOOLS COMPUTATION ON JUNCTIONS
# Then to compute the mean read depth
#dup_save$MeanDepth1 <- (dup_save$ReadDepthPool1/dup_save$SizePool1*150)*2
#dup_save$MeanDepth2 <- (dup_save$ReadDepthPool2/dup_save$SizePool1*150)*2

# Compute the frequencies
#dup_save$freqLL <- dup_save$NormReadsPool1/dup_save$MeanDepth1
#dup_save$freqLF <- dup_save$NormReadsPool2/dup_save$MeanDepth2
#dup_save_bis <- dup_save

# And finally, using the frequencies, compute the FST for each set of allelic frequencies
dup_save_bis$var <- 1/2*(dup_save_bis$freqLL^2+dup_save_bis$freqLF^2) -
  ((dup_save_bis$freqLL+dup_save_bis$freqLF)/2)^2
dup_save_bis$Fst <- dup_save_bis$var /
  ((dup_save_bis$freqLL+dup_save_bis$freqLF)/2 * (1-(dup_save_bis$freqLL+dup_save_bis$freqLF)/2))

hist(dup_save_bis$Fst, breaks=100)
max(dup_save_bis$Fst, na.rm=T) # 5.517041
View(dup_save_bis)

summary(as.factor(dup_save_bis$Type)) # only dupB if we do DiffSupIns threshold
summary(as.factor(dup_save_bis$TypeDup)) # only dupB if we do DiffSupIns threshold
summary(as.factor(dup_save_bis$extra)) # only dupB if we do DiffSupIns threshold

library(xlsx)
write.xlsx(dup_save_bis, file="dup_save_withcov.xls", 
           row.names = F)
```

After computing the Fst, compute informations about frequencies:

```{r}
library(xlsx)
dup_save_bis <- read.xlsx(file="dup_save_withcov.xls", sheetName = "Sheet1")

# Create a column with the information about frequency
dup_save_bis$fq_info <- NA

dup_save_bis$fq_info[dup_save_bis$freqLL==0 
                     & dup_save_bis$Type=="everted"] <- "spe_LF_duptand" 
# avec dup_save_bis$extra=="DUP_INV" : 556
dup_save_bis$fq_info[dup_save_bis$freqLF==0 
             & dup_save_bis$Type=="everted"] <- "spe_LL_duptand" # 3
dup_save_bis$fq_info[dup_save_bis$freqLF!=0 & dup_save_bis$freqLL!=0
             & dup_save_bis$freqLF>dup_save_bis$freqLL
             & dup_save_bis$Type=="everted"] <- "LF_sup_duptand" # 1306
dup_save_bis$fq_info[dup_save_bis$freqLF!=0 & dup_save_bis$freqLL!=0
             & dup_save_bis$freqLF<dup_save_bis$freqLL
             & dup_save_bis$Type=="everted"] <- "LL_sup_duptand" # 204
dup_save_bis$fq_info[dup_save_bis$freqLF!=0 & dup_save_bis$freqLL!=0
             & dup_save_bis$freqLF==dup_save_bis$freqLL
             & dup_save_bis$Type=="everted"] <- "same_fq_duptand" # 2

dim(dup_save_bis[dup_save_bis$Type=="everted",]) # 3725

dup_save_bis$fq_info[dup_save_bis$freqLL==0 
                     & dup_save_bis$Type=="inverted"] <- "spe_LF_dupinv" # 685
dup_save_bis$fq_info[dup_save_bis$freqLF==0 
                     & dup_save_bis$Type=="inverted"] <- "spe_LL_dupinv" # 4
dup_save_bis$fq_info[dup_save_bis$freqLF!=0 & dup_save_bis$freqLL!=0
                     & dup_save_bis$freqLF>dup_save_bis$freqLL
                     & dup_save_bis$Type=="inverted"] <- "LF_sup_dupinv" # 1936
dup_save_bis$fq_info[dup_save_bis$freqLF!=0 & dup_save_bis$freqLL!=0
                     & dup_save_bis$freqLF<dup_save_bis$freqLL
                     & dup_save_bis$Type=="inverted"] <- "LL_sup_dupinv" # 255
dup_save_bis$fq_info[dup_save_bis$freqLF!=0 & dup_save_bis$freqLL!=0
                     & dup_save_bis$freqLF==dup_save_bis$freqLL
                     & dup_save_bis$Type=="inverted"] <- "same_fq_dupinv" # 4

dim(dup_save_bis[dup_save_bis$Type=="inverted",]) # 5122

# Create a column with the information about frequency

dup_save_bis$readratio_info <- NA

dup_save_bis$readratio_info[dup_save_bis$ReadDepthRatio > dup_save_bis$LowerThreshold 
                         | dup_save_bis$ReadDepthRatio < dup_save_bis$UpperThreshold] <- "between"
dup_save_bis$readratio_info[dup_save_bis$ReadDepthRatio > dup_save_bis$UpperThreshold] <- "up"
dup_save_bis$readratio_info[dup_save_bis$ReadDepthRatio < dup_save_bis$LowerThreshold] <- "low"

summary(as.factor(dup_save_bis$readratio_info)) 
# between  low  up 
# 7112     986  749

p <- ggplot(dup_save_bis, aes(as.factor(fq_info), fill=readratio_info)) +
  geom_bar(stat="count", color="black", position=position_dodge())+
  geom_text(aes(label=..count..), stat="count", vjust=0)
print(p)
# Graphique fini
p + labs(x="Frequency information", y = "Number of events")+
  theme_classic()

# Probleme : "FP" labelled si cov >0 !!! pas normalisé...
# We then can try with other thresholds: q25/75
q25 <- read.xlsx(file="quantiles_25.xlsx", sheetName = "Feuil1")
dup_save_bis <- merge(dup_save_bis, q25, by="SizeClass")
dup_save_bis$q25 <- as.numeric(as.character(dup_save_bis$q25))
dup_save_bis$q75 <- as.numeric(as.character(dup_save_bis$q75))

inv_bis <- dup_save_bis[dup_save_bis$Type=="inverted" &
                          dup_save_bis$freqLF>dup_save_bis$freqLL & 
                          dup_save_bis$ReadDepthRatio>dup_save_bis$q25 &
                          dup_save_bis$ReadDepthRatio<dup_save_bis$q75,] # 1 405
inv_bis$TypeEvent <- "inversion"

dup_inv_bis <- dup_save_bis[dup_save_bis$Type=="inverted" &
                              dup_save_bis$freqLF>dup_save_bis$freqLL & 
                              dup_save_bis$ReadDepthRatio<dup_save_bis$q25,] # 692
dup_inv_bis$TypeEvent <- "inverted_dup"

dup_tand_bis <- dup_save_bis[dup_save_bis$Type=="everted" &
                               dup_save_bis$freqLF>dup_save_bis$freqLL &
                               dup_save_bis$ReadDepthRatio<dup_save_bis$q25,] # 538
dup_tand_bis$TypeEvent <- "tandem_dup"

quantile(inv_bis$Fst, na.rm=T, probs = 0.90) # 90% = 0.0396994
quantile(dup_inv_bis$Fst, na.rm=T, probs = 0.90) # 90% = 0.04014379
quantile(dup_tand_bis$Fst, na.rm=T, probs = 0.90) # 90% = 0.03636364

# compute empirical pval en FST
n=dim(inv_bis)[1]
library(dplyr)
inv_bis = inv_bis %>% arrange(desc(Fst)) 
inv_bis$Fst_pval=rep(1,n)
# calculer cb de valeurs sont plus petites que celle la
inv_bis$Fst_pval = (1:n)/n 

n=dim(dup_inv_bis)[1]
dup_inv_bis = dup_inv_bis %>% arrange(desc(Fst)) 
dup_inv_bis$Fst_pval=rep(1,n)
# calculer cb de valeurs sont plus petites que celle la
dup_inv_bis$Fst_pval = (1:n)/n 
View(dup_inv_bis)

all_clean <- rbind(dup_inv_bis,dup_tand_bis)
all_clean <- rbind(all_clean,inv_bis)
all_clean$fq_info <- substring(all_clean$fq_info, 1,6)
all_clean <- all_clean[all_clean$scaffold!="NA",]

# Add linkage groups
library(xlsx)
LG <- read.xlsx(file="scaff_LG_170122.xls", sheetName = "Sheet1",
                header = T)
LG$scaffold <- paste0(LG$seqid, ".1")

all_clean <- merge(all_clean, LG[,c(4,6)], 
                 by=c("scaffold"), all.x = T)

library(ggplot2)
p <- ggplot(all_clean, aes(as.factor(fq_info), fill=TypeEvent)) +
  geom_bar(stat="count", color="black", position=position_dodge())+
  geom_text(aes(label=..count..), stat="count", vjust=0)
print(p)
# Graphique fini
p+labs(x="Frequency information", y = "Number of events")+
  theme_classic()
```

Are reads delineating those events falling inside TE positions ? So could it be false positives events ?

```{r}
blast_TE=read.csv("results-blast-TE.txt", header=F, sep="\t") # 466,713
colnames(blast_TE) <- c("id_query","id_target","sq_identity","alignment_length",
                        "nb_mismatches","nb_gap_openings","domain_start_query",
                        "domain_end_query","domain_start_target",
                        "domain_end_target","evalue","bit_score")
blast_TE <- blast_TE[blast_TE$sq_identity > 0.8,] # 379,736
head(blast_TE) # we need id_target, domain_start_target, domain_end_target
blast_TE <- blast_TE[,c(2,9,10)]

# check in this table: how many events from scaffolds enriched in SVs AND in TE ?
# with at least 100 bp of overlap on intervals

library(GenomicRanges)
library(IRanges)
ranges_TE <- GRanges(seqnames=blast_TE$id_target, 
                     ranges=IRanges(blast_TE$domain_start_target, blast_TE$domain_end_target))
ranges_SV_start <-GRanges(seqnames=all_clean$scaffold, 
                     ranges=IRanges(all_clean$start_event, all_clean$start_event_bis))
ranges_SV_end <-GRanges(seqnames=all_clean$scaffold, 
                        ranges=IRanges(all_clean$end_event_bis, all_clean$end_event))

ranges_common_start <- subsetByOverlaps(ranges_SV_start, ranges_TE, minoverlap=100) # 393
length(ranges_common_start)
ranges_common_end <- subsetByOverlaps(ranges_SV_end, ranges_TE, minoverlap=100) # 414
length(ranges_common_end)

ranges_SV_start <- as.data.frame(ranges_SV_start)
ranges_common_start <- as.data.frame(ranges_common_start)
ranges_SV_start_sub <- all_clean[paste0(all_clean$scaffold, all_clean$start_event) 
                                 %in% paste0(ranges_common_start$seqnames, ranges_common_start$start),] # 393

ranges_SV_end <- as.data.frame(ranges_SV_end)
ranges_common_end <- as.data.frame(ranges_common_end)
ranges_SV_end_sub <- all_clean[paste0(all_clean$scaffold, all_clean$end_event_bis) 
                               %in% paste0(ranges_common_end$seqnames, ranges_common_end$start),] # 414

ranges_common_both <- ranges_SV_start_sub[ranges_SV_start_sub$Coords %in% ranges_SV_end_sub$Coords,] # 138
dim(ranges_common_both) # 138
View(ranges_common_both)

all_clean <- all_clean[!(all_clean$Coords %in% ranges_common_both$Coords),] # 2 497

library(plyr)
count(as.factor(all_clean$TypeEvent))
# inversion	= 1338			
# inverted_dup	= 650			
# tandem_dup = 509	

write.xlsx(all_clean, file="all_clean.xls", 
           row.names = F)

all_clean$freqLF[all_clean$freqLF>1] <- 1
all_clean$freqLL[all_clean$freqLL>1] <- 1
```


```{r}
# Find whether some events fall within genes
library(xlsx)
merge_gff <- read.xlsx(file="merge_gff.xlsx", sheetName = "Sheet1")
annot_R <- read.xlsx2(file="annotation_bedbug_221121.xls",
                      sheetName = "insecticide_resistance_100122")
merge_gff <- merge(merge_gff, annot_R[,c(1,7)], by="Name", all.x=T)

all_clean <- read.xlsx(file="all_clean.xls", sheetName = "Sheet1") 
test <- all_clean[all_clean$scaffold=="NW_019392721.1" |
                    all_clean$scaffold=="NW_019392763.1" |
                    all_clean$scaffold=="NW_019942502.1",]
View(test)

final_dup <- merge(all_clean, merge_gff, by.x="scaffold", by.y="seqid", all.x = T) # 142 480
colnames(final_dup)[46] <- "start_gene"
colnames(final_dup)[47] <- "end_gene"

final_dup$start_event <- as.numeric(as.character(final_dup$start_event))
final_dup$end_event <- as.numeric(as.character(final_dup$end_event))

# Then, we want to select lines where genes and events are overlapping

sub_dup <- final_dup[(final_dup$start_event-1000 < final_dup$start_gene &
                       final_dup$end_event+1000 > final_dup$end_gene),] # 1kb around events: 7 457 events in genes 
sub_dup <- sub_dup[!is.na(sub_dup$TypeDup),] # 7 426
length(unique(sub_dup$Name)) # 4118 genes

summary(as.factor(sub_dup$fq_info))
# sub_dup <- sub_dup[,c(1,9,11,33:36,38,43:49)] 

library(xlsx)
write.xlsx(sub_dup, file="sub_dup.xls", row.names = F) 

### Load now
sub_dup <- read.xlsx(file="sub_dup.xls", sheetName = "Sheet1") 

length(unique(sub_dup$Name)) # 4118 genes BUT a single gene could be in different TypeEvent

sub_dup_outliers <- sub_dup[((sub_dup$TypeEvent=="inversion" & sub_dup$Fst>=0.0396994) |
          (sub_dup$TypeEvent=="inverted_dup" & sub_dup$Fst>0.04014379) |
          (sub_dup$TypeEvent=="tandem_dup" & sub_dup$Fst>0.03636364)),]
```

Then, we want to extract only resistance genes:

```{r}
resist_dup <- sub_dup[!is.na(sub_dup$category),] # 255
summary(as.factor(resist_dup$TypeEvent))
# inversion inverted_dup   tandem_dup 
#       134           86           35 
         
p <- ggplot(resist_dup, aes(as.factor(fq_info), fill=TypeEvent)) +
  geom_bar(stat="count", color="black", position=position_dodge())+
  geom_text(aes(label=..count..), stat="count", vjust=0)
print(p)
p+labs(x="Frequency information", y = "Number of events")+
  theme_classic()

dim(resist_dup)[1] # 255
length(unique(resist_dup$Coords)) # 143 events in resistance genes
length(unique(resist_dup$Name)) # 130 genes

### SAVE ###
library(xlsx)
write.xlsx(resist_dup, file="resist_dup.xls", row.names = F) 
############

resist_dup <- read.xlsx(file="resist_dup.xls", sheetName = "Sheet1")
View(resist_dup) 

# We want to subset resist_dup according to scaffolds overexpressed
# NW_019392665.1 dup_inverted or dup_tandem or inversion
# NW_019392673.1 inversion
# NW_019392686.1 inversion
# NW_019392721.1 inversion or dup_inverted
# NW_019392754.1 dup_inverted or dup_tandem
# NW_019392785.1 inversion
# NW_019392813.1 inversion

resist_dup_scaffover <- resist_dup[
  resist_dup$scaffold=="NW_019392665.1" & (resist_dup$TypeEvent=="inversion" | resist_dup$TypeEvent=="inverted_dup" | resist_dup$TypeEvent=="tandem_dup") |
  resist_dup$scaffold=="NW_019392673.1" & resist_dup$TypeEvent=="inversion" |
  resist_dup$scaffold=="NW_019392686.1" & resist_dup$TypeEvent=="inversion" |
  resist_dup$scaffold=="NW_019392721.1" & (resist_dup$TypeEvent=="inversion"  | resist_dup$TypeEvent=="inverted_dup") |
  resist_dup$scaffold=="NW_019392754.1" & (resist_dup$TypeEvent=="inverted_dup" | resist_dup$TypeEvent=="tandem_dup") |
  resist_dup$scaffold=="NW_019392785.1" & resist_dup$TypeEvent=="inversion" |
  resist_dup$scaffold=="NW_019392813.1" & resist_dup$TypeEvent=="inversion",]
# 64 obs
View(resist_dup_scaffover)

# Then, we want to know how many of each event type fall in R genes
resist_dup_scaffover <- resist_dup_scaffover[,c(1,43:45,48)]

resist_dup_scaffover <- resist_dup_scaffover %>% 
  group_by(scaffold,TypeEvent,LG,Name,product) %>% 
  dplyr::summarize(count = n())

library(dplyr)
length(unique(resist_dup_scaffover$Name)) # 16 genes
resist_dup_scaffover <- as.data.frame(resist_dup_scaffover)

library(xlsx)
write.xlsx(resist_dup_scaffover, file="scaffxtrem_resist.xls")

resist_inv <- resist_dup[resist_dup$TypeEvent=="inversion",] # 134
length(unique(resist_inv$Name)) # 80 events

resist_dup_inv <- resist_dup[resist_dup$TypeEvent=="inverted_dup",] # 86
length(unique(resist_dup_inv$Name)) # 50 events

resist_dup_tand <-resist_dup[resist_dup$TypeEvent=="tandem_dup",] # 24
length(unique(resist_dup_tand$Name)) # 30 events

# Then filter on FST
resist_dup_tand <- resist_dup_tand[resist_dup_tand$Fst>=0.03636364,] # 0 obs
resist_dup_inv <- resist_dup_inv[resist_dup_inv$Fst>=0.04014379,] # 9 obs
resist_inv <- resist_inv[resist_inv$Fst>=0.0396994,] # 10 obs

outliers_resist <- rbind(resist_dup_inv, resist_dup_tand)
outliers_resist <- rbind(outliers_resist, resist_inv)

length(unique(resist_dup$Coords)) # 143 events
length(unique(resist_dup$Name)) # 130 genes

write.xlsx(outliers_resist, file="outliers_resist.xls")
```




