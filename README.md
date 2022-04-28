![Image 1](bedbugs.png)

Contact : chloe.haberkorn@univ-lyon1.fr / chloehbk@gmail.com

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
	- [Alternative alleles](#Alternative-alleles)
	- [Synonymous or not](#Synonymous-or-not)

 - **[Copy Number Variation](#Copy-Number-Variation)**
 	- [Computing average coverage by gene](#Computing-average-coverage-by-gene)
	- [Selecting amplified genes](#Selecting-amplified-genes)



## Installing tools

Here are the tools and versions used (on a cluster): 
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

Raw sequences (fastq.gz files) are available on SRA: xx
```
mkdir /your-path/PoolSeq_Clec/Raw_Clec
```

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

We can compute coverage per base (i.e. number of reads mapping at a position of the reference genome)

```
/your-path/Tools/Tools/bedtools2/bin/bedtools genomecov -ibam SF_mapped_sorted.bam -d > SF_cov.txt
```

We choose to exclude of our analysis coverage over 50 bp, which corresponds to >95% quantile for all populations, in order to avoid bias due to very high coverage. Coverages computed here were also used for [Copy Number Variation](#Copy-Number-Variation) analysis.



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

data <- data[!(data$contig=="NC_030043.1"),] # Remove mitochondrial : 4,251,923 SNPs
``` 
The final dataset contains 4 251 923 loci. However, depending on the comparison the number of SNP conserved differ because of the maf filter (applied for each pairwise comparison).


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
# 2,175,417 SNPs

true_ref <- (data$ref == data$nucleotides) & (data$alt !=  data$nucleotides)
# 2,075,247 SNPs

# Check whether there are alleles of ref & alt != alleles of the reference genome (harlan strain)
wrong_all <- (data$ref != data$nucleotides) & (data$alt !=  data$nucleotides)
# yes indeed ! 1,259 SNPs

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

### Alternative alleles


On veut identifier les SNPs où LL a l’allele de ref Harlan, et LF l’allele alt

```{r}
NT <- read.table("~/Desktop/data_clean_nt.txt", sep="\t", header=TRUE)

head(NT)
# on veut LL ref > 1/2 reads (>1 indiv sur 2 chez LL porte l'allele de ref)
# & LF ref < 1/2 reads (>1 indiv sur 2 chez LF porte l'allele alt)
test <-  NT[(NT$LL_ref > NT$LL_tot/2) & (NT$LF_ref < NT$LF_tot/2),] # 298 788 obs
298788/5574186

test <-  NT[(NT$LL_ref == 1) & (NT$LF_ref == 0),] # 20 491 obs
```

### Synonymous or not

## Copy Number Variation

### Computing average coverage by gene

We computed coverage per base using bedtools, as in [Analysing coverage](#Analysing-coverage).

Using the gff, we created a table to keep only exonic part, with gene annotation:
```
gtffile <- file.path("/your-path/GCF_000648675.2_Clec_2.1_genomic.gff") 
library(ape)
gff=read.gff(gtffile)

scaff_exons=subset(gff,type == "exon") # 28 5072 rows
scaff_exons=subset(scaff_exons,source != "RefSeq") # 28 5047 rows

library(tidyverse)
scaff_exons = separate(data=scaff_exons, 
                      into=c("ID","gene"),
                      col=attributes, sep="gene=")

scaff_exons = separate(data=scaff_exons, 
                       into=c("gene","drop"),
                       col=gene, sep=";",extra="merge")

scaff_exons = separate(data=scaff_exons, 
                      into=c("drop2","product"),
                      col=drop, sep="product=", extra="merge")

scaff_exons = separate(data=scaff_exons, 
                       into=c("annotation","drop3"),
                       col=product, sep=";", extra="merge")

scaff_exons <- scaff_exons %>%
  select(seqid, gene, annotation)
scaff_exons <- unique(scaff_exons) # 27 760 rows

library(xlsx)
write.xlsx2(scaff_exons, file="/your-path/scaff_exons.xls", sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)

# Remove manually all "transcript variant X..." and "%2C" and create a second sheet ("Feuil1") without Trna

scaff_exons = read.xlsx("/your-path/scaff_exons.xls", sheetName = "Feuil1", header=T)
scaff_exons <- unique(scaff_exons) # 10 494 rows 
```

We compute a file with a position per base for each exon:
```
/your-path/Tools/R-3.5.2/bin/R

scaff_exon=read.table("scaff_exon.txt",header=T)

library(dplyr)
library(tidyr)

exon_perbase <- 
  scaff_exon %>% 
  rowwise %>% 
  mutate(pos = paste(seq(start, 
                         end), 
                     collapse = ",")) %>%
  tidyr::separate_rows(pos) %>%
  mutate(pos = as.integer(pos))
```

Coverages were merged by exon for each strain (example here with London Field):
```
library(dplyr)
library(tidyr)
map_LF=read.table("LF_cov.txt") 
colnames(map_LF) <- c("seqid","pos","cov")

genes_cov <- right_join(exon_perbase, map_LF)
# > Joining, by = c("seqid", "pos")

exons_cov_LF <- genes_cov %>% 
  group_by(seqid, start, end, gene) %>%
  summarise(sum_coverage = sum(cov))
```

We merged exons with the same Gene ID together to build genes. To do so, cumulated length of exons inside a gene was computed, as well as cumulated coverage (example here with London Field).
```
exons_cov_LF$length_exon <- exons_cov_LF$end - exons_cov_LF_clean$start

genes_cov_LF <- exons_cov_LF %>% 
  group_by(seqid, gene) %>% 
  summarise_at(c("length_exon","sum_coverage"),sum)
colnames(genes_cov_LF) <- c("seqid", "gene", "length_gene", "cov_gene")

```

Now we can add annotation:
```
genes_cov_LF <- merge(genes_cov_LF, scaff_exons, by=c("gene","seqid"), # using both gene and seqid because Trna had same name on different scaffold
                      all.x = F, all.y = F) # 27 824 for each
```

Keep only one product column (to deal with "protein Malvolio-like" for example):
```
genes_cov_LF_clean <- genes_cov_LF %>% 
  group_by(seqid, gene, length_gene) %>% 
  filter(annotation==min(annotation) | is.na(annotation)) # 9 947 tows

```

In "sum_coverage" we have the "cumulated coverage" along the gene for each gene.
Which gave us: average coverage by base pair = sum_coverage/length of gene
```
genes_cov_LF_clean$avg_coverage <- genes_cov_LF_clean$cov_gene/genes_cov_LF_clean$length_gene
```


### Selecting amplified genes

To detect amplified genes, we combined a density-based approach with a linear model and a ratio of average coverage between London Lab and London Field.

We first plot those average coverages using geom density 2D (using tutorials available here: https://r-charts.com/correlation/contour-plot-ggplot2/ and http://slowkow.com/notes/ggplot2-color-by-density/).
```
sub_cov_genes <- genes_cov_LL_clean[c(6,7,9)] # keep only LL and LF avg_coverage and XXXX

library(ggplot2)
gg <- ggplot(data = sub_cov_genes, aes(x = avg_coverage, y = avg_coverage_LF)) +
  geom_point(cex = 0.5, col="grey") +
  xlim(0,100) + ylim(0,100) +
  geom_density_2d_filled(alpha = 0.5, bins=15)+
  theme_bw()
gg 
```

Levels denoted by contours were extracted by going into the  ggplot build object, and a layer that relies on given contour level was added:
```
gb <- ggplot_build(gg)
contour_levels <- unique(gb[["data"]][[2]][["level"]])

get_density <- function(x, y, n = 9947) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

sub_cov_genes$density <- get_density(sub_cov_genes$avg_coverage,
                                     sub_cov_genes$avg_coverage_LF)
summary(sub_cov_genes$density<0.0008) # 2599 values
```

Data were trimmed according to quantile 1% and 99%:
```
quantile(genes_cov_LL_clean$avg_coverage, probs = c(0.01,0.99)) 
quantile(genes_cov_LF_clean$avg_coverage, probs = c(0.01,0.99)) 

sub_cov_genes_clean <- sub_cov_genes[sub_cov_genes$avg_coverage > 11.48672 &
                        sub_cov_genes$avg_coverage < 126.76237 &
                        sub_cov_genes$avg_coverage_LF > 13.22879 &
                        sub_cov_genes$avg_coverage_LF < 148.84460,] 
```

We compute a linear regression between average coverage of London Field and London Lab:
```
lm_cov <- lm(sub_cov_genes_clean$avg_coverage_LF ~ sub_cov_genes_clean$avg_coverage)
```

To highlight points above ratio 1.5 for this linear regression, we first added a column to our data, and then plot it:
```
sub_cov_genes_clean$ratio_lm <- 1.5*(lm_cov$coefficients[2]*sub_cov_genes_clean$avg_coverage)+lm_cov$coefficients[1]

library(viridisLite)
gg2 <- ggplot(data = sub_cov_genes_clean, aes(x = avg_coverage, y = avg_coverage_LF)) +
  geom_point(aes(avg_coverage, avg_coverage_LF, color = density), 
             cex=0.2) + 
  scale_color_viridis() +
  geom_point(data = sub_cov_genes_clean %>% filter(density <= 0.0008), color = "grey", size = 0.2) +
  geom_point(data = sub_cov_genes_clean %>% 
               filter(density <= 0.0008 & avg_coverage_LF > ratio_lm), 
             color = "black", size = 0.2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black", size=0.3)) +
  geom_vline(xintercept=11.48672, size = 0.2, linetype="dashed") +
  #geom_vline(xintercept=126.76237, size = 0.2, linetype="dashed") +
  geom_hline(yintercept=13.22879, size = 0.2, linetype="dashed") +
  #geom_hline(yintercept=148.84460, size = 0.2, linetype="dashed") +
  geom_abline(intercept = lm_cov$coefficients[1] , slope = lm_cov$coefficients[2]*1.5, size = 0.3, col="blue4")+
  xlim(11.48672,140) + 
  ylim(13.22879,140)
```

Finally, we extracted those 25 points which correspond to our final set of 25 CNV, and added all informations:
```
set_cnv <- sub_cov_genes_clean[sub_cov_genes_clean$avg_coverage_LF > sub_cov_genes_clean$ratio_lm & 
                              sub_cov_genes_clean$density <= 0.0008,] # 25 values
			      
set_cnv <- merge(set_cnv,genes_cov_LL_clean, all.x = T, all.y = F, by=c("avg_coverage","avg_coverage_LF","ratio"))

LG <- read.xlsx(file="/your-path/scaff_LG.xls", sheetName = "Sheet1",
                header = T)
common <- intersect(LG$seqid, set_cnv$seqid) # 12 values (13, but two on NW_019392676.1)

set_cnv <- merge(set_cnv, LG, by=c("seqid"), all.x = T)

library(xlsx)
write.xlsx(set_cnv, file="/your-path/set_cnv.xls", 
           row.names = F)
```

