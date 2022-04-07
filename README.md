![Image 1](bedbugs.png)

Contact : chloe.haberkorn@univ-lyon1.fr

### Table of Contents

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
	- [Alternative alleles](#Alternative-alleles)




## Installing tools

Here are the tools and versions used (on a cluster): 
- FastQC 
- Trimmomatic v0.39
- FastUniq v1.1
- BWA v0.74
- Samtools v1.9
- Bedtools v2.29.1
- PoPoolation 2 v1201
- BayPass v2.2
- R v3.5.2

They will be store in /your-path/Tools.
We also used R on a computer with packages poolfstat v2.0.0, pcadapt v4.3.3, VariantAnnotation v1.34.0, and GenomicFeatures v1.40.1.




## Pool-seq data processing

The goal is first to map *Cimex lectularius* PoolSeq samples (London Lab, London Field, German Lab and Sweden Field - pools of 30 individuals) on reference genome.
We used the recent reference genome and annotation, avalaible here: https://www.ncbi.nlm.nih.gov/assembly/GCF_000648675.2

We will have to download a few softs.

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

We can compute coverage per base (i.e. number of reads mapping at a position of the reference genome)

```
/your-path/Tools/Tools/bedtools2/bin/bedtools genomecov -ibam SF_mapped_sorted.bam -d > SF_cov.txt
```

We choose to exclude of our analysis coverage over 50 bp, which corresponds to >95% quantile for all populations, in order to avoid bias due to very high coverage. Coverages computed here were also used for copy number variation analysis.




## Overall SNPs analyzes

The goal was to understand what differenciates the four *Cimex lectularius* PoolSeq samples: London Lab, London Field, German Lab, Sweden Field. 
Our hypothesis was that we could be able to find candidate loci correlated with their insecticide resistance phenotypes - resistant for Field strains and susceptible for Lab strains. 

For following analysis, we excluded the scaffold "NC_030043.1", which corresponds to the mitochondrial genome. Indeed, for one copy of the nuclear genome, there are several copies of the nuclear genome. Furthermore, the mitochondrial genome does not evolve like the nuclear genome (not the same mutation rate, no recombination, maternal transmission). 

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
                                min.maf = 0.1)
# pool.index = indexes of the pools (at least two), that should be selected to create the new pooldata objec
# min.maf correspond to minimal allowed Minor Allele Frequency
# max.cov of 50 : corresponds to coverage >q95
# Data consists of 4,251,927 (4,25 M) SNPs for 4 Pools
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

### Computing FST

For each SNP, we can compute SNP-specific pairwise FST for each comparisons between strains (GL_vs_LF, GL_vs_SF, GL_vs_LL, LF_vs_SF, LF_vs_LL, SF_vs_LL), thanks to the option "output.snp.values = TRUE": 
``` 
PairwiseFST_all = na.omit(computePairwiseFSTmatrix(pooldata_sub, method = "Anova",
						   min.cov.per.pool = -1, 
                                                   max.cov.per.pool = 1e+06,
                                                   min.maf = -1, 
                                                   output.snp.values = TRUE))
# Warning : here, min.maf by pairs, not for the 4 strains together !
``` 

Then, we merged informations together for following analysis:
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

data <- data[!(data$contig=="NC_030043.1"),] # Remove mitochondrial : 4,251,923 SNPs
``` 
The final dataset contains 4 251 923 loci. However, depending on the comparison the number of SNP conserved differ because of the maf filter (applied for each pairwise comparison).

We built a table combining all informations:
``` 
baypass=read.table("/Users/chloe/Documents/Cluster/baypass_220321_results.txt", sep=" ")
colnames(baypass) <- c("contig","position","C2_std","LOG_C2","BF.dB","XtXst","LOG_xtx")

data <- merge(data, baypass, by.x=c("contig","position"), by.y=c("SCAFFOLD","POSITION"), all.y=F) # add baypass to SNP/FST informations
``` 

Since reference allele is "chosen arbitrarily in each pool" by poolfstat, we need to correct it. Indeed, we want the reference allele to match the reference genome allele, i.e. Harlan strain allele.

#1 Add a column reference allele related to Harlan.

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

mv data.fa.out data_out.bed
```
! Check that there is no number in scientific writing in BED entry ! 

```{r}
# Verifier que l'output de bedtools correspond à ce qu'on observait dans le tableau data_all

data_out <- read.table("~/Desktop/data_out.bed", sep="\t")
dim(data_out)

chr <- seq(from = 1, to = 11148372, by = 2)
nt <- seq(from = 2, to = 11148372, by = 2)

library(dplyr)
positions <- filter(data_out, row_number() %in% chr)
nucleotides <- filter(data_out, row_number() %in% nt)

tab_out <- data.frame(positions, nucleotides)
colnames(tab_out) <- c("positions","nucleotides")

library(tidyr)
tab_out <- transform(tab_out, pos=do.call(rbind, strsplit(positions, ':', fixed=TRUE)), stringsAsFactors=F)
tab_out <- transform(tab_out, position=do.call(rbind, strsplit(pos.2, '-', fixed=TRUE)), stringsAsFactors=F)
tab_out <- transform(tab_out, scaff=do.call(rbind, strsplit(pos.1, '>', fixed=TRUE)), stringsAsFactors=F)

tab_out_ok <- tab_out[,c("nucleotides","position.2","scaff.2")]
head(data_bed) # checkcolumns names
colnames(tab_out_ok) <- c("nucleotides","position","contig")
head(tab_out_ok)
tab_out_ok$position <- as.numeric(as.character(tab_out_ok$position))

# merge tab_out_ok et data_all
data_all$position <- as.numeric(as.character(data_all$position))
data_all = data_all %>% arrange(position)

head(data_all)
head(tab_out_ok) # verifier que les tableaux sont ordonnes pareil avant le cbind

data_nt <- cbind(data_all, tab_out_ok$nucleotides)
head(data_nt)
colnames(data_nt)[52] <- "nucleotides"
```


Pour obtenir les alleles de ref :

```{r}
head(data_nt) # nucleotide de ref et alt en majuscule

library(stringr)
# mettre tous les nucleotides en majuscule
data_nt$nucleotides <- str_to_upper(data_nt$nucleotides) 

# as numeric
data_nt$LL_ref <- as.numeric(as.character(data_nt$LL_ref))
data_nt$LL_tot <- as.numeric(as.character(data_nt$LL_tot))
data_nt$LF_ref <- as.numeric(as.character(data_nt$LF_ref))
data_nt$LF_tot <- as.numeric(as.character(data_nt$LF_tot))
data_nt$GL_ref <- as.numeric(as.character(data_nt$GL_ref))
data_nt$GL_tot <- as.numeric(as.character(data_nt$GL_tot))
data_nt$SF_ref <- as.numeric(as.character(data_nt$SF_ref))
data_nt$SF_tot <- as.numeric(as.character(data_nt$SF_tot))

NT = as.data.table(data_nt)
# si difference allele "ref" =/= "nucleotide du genome de ref" et allele "alt" == "nucleotide"
wrong_ref <- (NT$ref != NT$nucleotides) & (NT$alt ==  NT$nucleotides) 

wrong_ref_count <- NT[(NT$ref != NT$nucleotides) & (NT$alt ==  NT$nucleotides),] 
# 2 877 733 alleles "ref" qui sont en realite des alleles alternatifs





# problème : certaines lignes ont un snp_count de zéro alors qu'elles sont justement autour d'un SNP ?
# explication : snp_count >2 pour calcul pi dans popoolation

snp_zero <- NT[(NT$snp_count_LL == 0) & (NT$snp_count_LF ==  0) & (NT$snp_count_GL ==  0) & (NT$snp_count_SF ==  0),] 
# 18311 
View(snp_zero)

snp_zero_all <- NT[(NT$snp_count_LL == 0) & (NT$snp_count_LF ==  0) 
                   & (NT$snp_count_GL ==  0) & (NT$snp_count_SF ==  0) 
                   & (NT$snp_count_LL_10kb == 0) & (NT$snp_count_LF_10kb ==  0) 
                   & (NT$snp_count_GL_w10kb ==  0) 
                   & (NT$snp_count_SF_w10kb ==  0),] # no value






true_ref_count <- NT[(NT$ref == NT$nucleotides) & (NT$alt !=  NT$nucleotides),] 
# 2 694 770  alleles "ref" qui sont les memes que l'allele du genome de ref

# verif si il y a des alleles de ref & alt != des allele du genome de ref (harlan strain)
wrong_all_count <- NT[(NT$ref != NT$nucleotides) & (NT$alt !=  NT$nucleotides),] 
# 1683 obs (wrong_ref_count + true_ref_count = 5572503 > 5574186 - 5572503 = 1683)
sum(is.na(NT$ref)) # clean !

# data_all_clean <- read.table("~/Desktop/data_all_clean.txt", sep="\t", header=TRUE)
# check if position is numeric before !
# data_all_clean = data_all_clean %>% arrange(position)
# data_nt = data_nt %>% arrange(position)
# head(data_all_clean)
# head(data_nt) # verifier que les tableaux sont ordonnes pareil avant le cbind
# data_clean_nt <- cbind(data_all_clean, data_nt$nucleotides)

# on remplace le count de ref par le count d'alt 
library(data.table)
NT <- NT[wrong_ref, LL_ref := LL_tot-LL_ref]
NT <- NT[wrong_ref, GL_ref := GL_tot-GL_ref]
NT <- NT[wrong_ref, LF_ref := LF_tot-LF_ref]
NT <- NT[wrong_ref, SF_ref := SF_tot-SF_ref]

# on remplace l'allele alt par l'allele ref et l'allele ref par le "vrai" allele ref
NT <- NT[wrong_ref, alt := ref]
NT <- NT[wrong_ref, ref := nucleotides]

head(NT)
head(data_nt) # compa pour verifier que ça marche
```


#2 Modifier aussi les fq  GL_p, LF_p, SF_p, LL_p

```{r}
# on recalcule les FA de l'allele de ref
NT$GL_p <- NT$GL_ref/NT$GL_tot
NT$LL_p <- NT$LL_ref/NT$LL_tot
NT$LF_p <- NT$LF_ref/NT$LF_tot
NT$SF_p <- NT$SF_ref/NT$SF_tot

write.table(NT, file = "~/Desktop/data_clean_nt.txt",quote=F, row.names=F, col.names=T, sep="\t")
```


#3 Ensuite, on veut identifier les SNPs où LL a l’allele de ref Harlan, et LF l’allele alt

```{r}
NT <- read.table("~/Desktop/data_clean_nt.txt", sep="\t", header=TRUE)

head(NT)
# on veut LL ref > 1/2 reads (>1 indiv sur 2 chez LL porte l'allele de ref)
# & LF ref < 1/2 reads (>1 indiv sur 2 chez LF porte l'allele alt)
test <-  NT[(NT$LL_ref > NT$LL_tot/2) & (NT$LF_ref < NT$LF_tot/2),] # 298 788 obs
298788/5574186

test <-  NT[(NT$LL_ref == 1) & (NT$LF_ref == 0),] # 20 491 obs

```



```{r}
set1=which(NT$LF_vs_LL > 0.4161941) # 52 502 avec q>99%
set3=which(NT$BF.dB>20 & NT$LOG_C2>3 & NT$LF_vs_LL > 0.4161941) # 65 avec q99
set4=which(test$BF.dB>20 & test$LOG_C2>3 & test$LF_vs_LL > 0.4161941) # 24 avec q99

set2=which(test$LF_vs_LL > 0.4161941) # 22 417 avec q>99%
candidates=test[unique(c(set2)),] 

setwd("~/Desktop")
write.table(candidates, file = "./Candidates_LLrefLFalt_99.txt",quote=F, row.names = F, col.names = T, sep="\t")

head(candidates)
```


## Selecting candidate SNPs

### Differentiated FST



### Selection with contrast between phenotypes

We performed a contrast analysis to identify SNPs associated with populations ecotypes. This trait (populations' ecotype) being binary, we can use C2 statistic (Olazcuaga et al., 2019) to identify those SNPs, rather than parametric models used to estimates Bayes' Factor (BF).

First, we converted Poolfstat SNPs data into BayPass input format:
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

We then merged two of output files together: poolfstatdata_220321_summary_contrast_snpdet.out and poolfstatdata_220321_summary_betai_reg_snpdet.out, in /your-path/PoolSeq_Clec/BayPass/baypass_220321_results.txt")

### Alternative alleles




