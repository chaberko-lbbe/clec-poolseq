![Image 1](bedbugs.png)

Contact : chloe.haberkorn@univ-lyon1.fr

### Table of Contents

- **[Mapping and processing Pool-seq genomes](#Mapping-and-processing-Pool-seq-genomes)**
	- [Install tools](#Install-tools)
	- [Getting the data](#Getting-the-data)
	- [Trimming](#Trimming)
	- [Removing duplicates](#Removing-duplicates)
	- [Mapping](#Mapping)
	- [Analysing coverage](#Analysing-coverage)

- **[Genetic differenciation of populations](#Genetic-differenciation-of-populations)**
	- [Install tools](#Install-tools)
	- [Detecting SNPs](#Detecting-SNPs)
	- [Compute FST](#Compute-FST)
	- [Estimate genetic polymorphism](#Estimate-genetic-polymorphism)

## Mapping and processing Pool-seq genomes

The goal is first to map *Cimex lectularius* PoolSeq samples on reference genome : London Lab, London Field, German Lab, Sweden Field (pools of 30 individuals).
We can benefit from the recent reference genome and annotation, avalaible here : https://www.ncbi.nlm.nih.gov/assembly/GCF_000648675.2

We will have to download a few softs.

### Install tools

BWA :
``` 
git clone https://github.com/lh3/bwa.git
cd bwa
make
/beegfs/data/chaberkorn/Tools/bwa/bwa index # Check that the soft is working 
```

Trimmomatic :
We will have to download the associated adapters. To do so : in fastq file, look at overrepresented sequences -> "TruSeq Adapter". These primers come from the TruSeq-3 library.
``` 
chaberkorn@pbil-deb:/beegfs/data/chaberkorn/Tools/trimmomatic-0.39
chaberkorn@pbil-deb:/beegfs/data/chaberkorn/Tools/TruSeq3-PE.fa
``` 

Bedtools :
``` 
chaberkorn@pbil-deb:/beegfs/data/chaberkorn/Tools/bedtools2
``` 

FastUniq :
``` 
chaberkorn@pbil-deb:/beegfs/data/chaberkorn/Tools/FastUniq
``` 

### Getting the data

```
#Create directory with raw sequences
mkdir /beegfs/data/chaberkorn/PoolSeq_Clec/Raw_Clec
cp -r /beegfs/data/varaldi/BEDBUGS/data_preliminary/*R*.fastq.gz /beegfs/data/chaberkorn/PoolSeq_Clec/Raw_Clec
tar xvf fastqc.tar.gz

#Create directory with reference genome 
mkdir /beegfs/data/chaberkorn/PoolSeq_Clec/Ref_Clec

#Open a terminal to move reference genome downloaded on computer to Beegfs
scp -r /Users/chloe/Documents/Cluster chaberkorn@pbil-deb:/beegfs/data/chaberkorn/PoolSeq_Clec/Ref_Clec
```

### Trimming

Keep default parameters :
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads 
  - TruSeq3-PE.fa = according to fastq file infos
  - 2 = seed mismatches
  - 30 = palindrome clip threshold
  - 10 = simple clip threshold
  - 2 = min Adapter Length
  - keepBothReads
LEADING:3
TRAILING:3 
SLIDINGWINDOW:4:20 -> too high, use 15 instead
MINLEN:50

```
cd /beegfs/data/chaberkorn/PoolSeq_Clec/Trimmed
nano trimm_LL.sh

#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem 10G
#SBATCH -t 24:00:00
#SBATCH -e /beegfs/data/chaberkorn/PoolSeq_Clec/Trimmed/Trimming_LL.error
#SBATCH -o /beegfs/data/chaberkorn/PoolSeq_Clec/Trimmed/Trimming_LL.out
#SBATCH -J Genome_trimming_Cimex_lectularius

DIR=/beegfs/data/chaberkorn/PoolSeq_Clec
DIRSOFT=/beegfs/data/chaberkorn/Tools
DIRFASTQ="$DIR"/Raw_Clec

/usr/local/jre1.8.0_202/bin/java -jar "DIRSOFT"/trimmomatic-0.39.jar PE -phred33 -trimlog LL_trim.log \
"DIRFASTQ"/LL_R1.fastq.gz "DIRFASTQ"/LL_R2.fastq.gz \
LL_R1_paired.fq.gz LL_R1_unpaired.fq.gz LL_R2_paired.fq.gz LL_R2_unpaired.fq.gz \
ILLUMINACLIP:"DIRSOFT"/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 

gunzip LF_*_paired.fq.gz

# Check filtering quality using fastq :
/beegfs/data/chaberkorn/Tools/FastQC/fastqc --java /usr/local/jre1.8.0_202/bin/java LF_R1_paired.fq
/beegfs/data/chaberkorn/Tools/FastQC/fastqc --java /usr/local/jre1.8.0_202/bin/java LF_R2_paired.fq
```

Count the number of sequences (starting with "@") with : 

```
grep -c @ LF_R2_paired.fq
```

### Removing duplicates

Using FastUniq

-t q : output sequence format - FASTQ format into TWO output files
-o : premier output
-p : deuxieme output
-c 1 : types of sequence descriptions for output - new serial numbers assigned by FastUniq (0 : the raw descriptions)

```
gunzip *_paired.fq.gz # R1 and R2 for each pop
# Create a text file with both input : R1 and R2 > input_LL.txt

nano LL_dup.txt

#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem 10G
#SBATCH -t 24:00:00
#SBATCH -e /beegfs/data/chaberkorn/PoolSeq_Clec/Trimmed/LL_dup.error
#SBATCH -o /beegfs/data/chaberkorn/PoolSeq_Clec/Trimmed/LL_dup.out
#SBATCH -J Genome_duplicates_Cimex_lectularius

/beegfs/data/chaberkorn/Tools/FastUniq/source/fastuniq -i input_LL.txt -t q -o LL_dup_1.fastq -p LL_dup_2.fastq -c 1

# Check how many duplicate sequences have been deleted:
/beegfs/data/chaberkorn/Tools/FastQC/fastqc --java /usr/local/jre1.8.0_202/bin/java LL_dup_1.fastq
/beegfs/data/chaberkorn/Tools/FastQC/fastqc --java /usr/local/jre1.8.0_202/bin/java LL_dup_2.fastq
```

### Mapping

```
nano mapping_LL.sh

#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem 10G
#SBATCH -t 24:00:00
#SBATCH -e /beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/mapping_LL.error
#SBATCH -o /beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/mapping_LL.out
#SBATCH -J Genome_mapping_Cimex_lectularius

/beegfs/data/chaberkorn/Tools/bwa/bwa mem -t 8 /beegfs/data/chaberkorn/PoolSeq_Clec/Ref_Clec/ncbi-genomes-2020-12-15/GCF_000648675.2_Clec_2.1_genomic.fna \
/beegfs/data/chaberkorn/PoolSeq_Clec/Trimmed/LL_dup_1.fastq /beegfs/data/chaberkorn/PoolSeq_Clec/Trimmed/LL_dup_2.fastq > LL_mapped.sam
```

Now we want to control mapping quality, while converting file from sam to bam

Two options : keep all reads (and know the percentage mapping via flagstat)
Or: separate mapped reads from unmapped reads (see /beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/SEP_MAP_UNMAP) :

```
/beegfs/data/soft/samtools-1.10/bin/samtools view -bS -F 4 GL_mapped.sam > GL_mapped.bam
/beegfs/data/soft/samtools-1.10/bin/samtools view -bS -f 4 SF_mapped.sam > SF_unmapped.bam
```

To keep only the reads whose mapping has a probability > 99% to be correct : add on samtools view "-q 20".

```
#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem 10G
#SBATCH -t 24:00:00
#SBATCH -e /beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/stats_SF.err
#SBATCH -o /beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/stats_SF.out
#SBATCH -J Genome_stats_Cimex_lectularius

/beegfs/data/soft/samtools-1.10/bin/samtools view -bS -@ 2 GL_mapped.sam > GL_mapped.bam  

/beegfs/data/soft/samtools-1.10/bin/samtools sort GL_mapped.bam -o GL_mapped_sorted.bam

rm GL_mapped.bam

/beegfs/data/soft/samtools-1.9/bin/samtools flagstat GL_mapped_sorted.bam > flagstat_GL.txt
```

### Analysing coverage 

Calculer les couvertures par base (combien de reads mappent à une position du genome de reference = 1 NT)
Pour uniquement les reads mappes (-F 4 ; mapped_sorted) 

```
nano cov_SF_mapped_sorted.sh

#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem 10G
#SBATCH -t 24:00:00
#SBATCH -e /beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/SEP_MAP_UNMAP/cov_SF_mapped_sorted.err
#SBATCH -o /beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/SEP_MAP_UNMAP/cov_SF_mapped_sorted.out
#SBATCH -J Genome_cov_Cimex_lectularius

/beegfs/data/chaberkorn/Tools/bedtools2/bin/bedtools genomecov -ibam SF_mapped_sorted.bam -d > SF_cov.txt
```

Then, we analyse them on R:

```
# Open R on the cluster: /beegfs/data/soft/R-3.5.2/bin/R
# Load the table :
map_SF=read.table("SF_cov_map.txt") 
head(map_SF)
summary(map_SF$V3)

# Create an additional column with the order of V1:
order <- seq(1:length(map_SF$V1))
map_SF$order <- order
head(map_SF)

map_SF_1 <- aggregate(V3~V1,data=map_SF,FUN=mean) # Aggregate with the mean of the Y column depending on X
map_SF_agg <- aggregate(order~V1, data=map_SF,FUN=mean) # Mean of column 'order'
map_SF <- map_SF_1
map_SF$order <- map_SF_agg$order # Add this column to data frame
map_SF <- map_SF[order(map_SF$order),] # Order the dataframe according to the averages of the 'order'
head(map_SF)

# How many scaffold with coverage=0 ?
cov_SF_null <- subset(map_SF,  map_SF$V3 == 0,)

write.table(map_SF, file="/beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/SEP_MAP_UNMAP/map_scaff_SF.txt")

map_GL_q=read.table("GL_cov_qual.txt")

# Try to find the basis differencialy covered between the two files :
test_GL <- setdiff(map_GL, map_GL_q) # find elements in x (the bigger file = map_GL) that are not in y (map_GL_q)

write.table(map_GL_q, file="/beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/SEP_MAP_UNMAP/map_scaff_GL_qual.txt")
```

Compute scaffold lengths:

```
map_LL=read.table("LL_cov_map.txt")
length_LL <- aggregate(V2~V1,FUN=length, data=map_LL)
write.table(length_LL, file="/beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/SEP_MAP_UNMAP/length_scaff_LL.txt")

# Merge tables and reorder by 'order' column:
map_LL=read.table("map_scaff_LL.txt")
length_LL=read.table("length_scaff_LL.txt")
map_length_LL <- merge(length_LL,map_LL, by.x="V1", by.y="V1")
colnames(map_length_LL) <- c("scaffold","length","mean_cov","order") 
map_length_LL <- map_length_LL[order(map_length_LL$order),]
write.table(map_length_LL, file="/beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/SEP_MAP_UNMAP/map_scaff_LL.txt")
```

Represent coverages on R:

```
setwd("~/Desktop/CloudStation/THESE/WholeGenome PoolSeq/Analyse coverage")

map_scaff_GL=read.table("map_scaff_GL.txt")
map_scaff_LL=read.table("map_scaff_LL.txt")
map_scaff_LF=read.table("map_scaff_LF.txt")
map_scaff_SF=read.table("map_scaff_SF.txt")

library(ggplot2)

map_scaff_GL$X <- 1:length(map_scaff_GL$scaffold)
map_scaff_LL$X <- 1:length(map_scaff_LL$scaffold)
map_scaff_LF$X <- 1:length(map_scaff_LF$scaffold)
map_scaff_SF$X <- 1:length(map_scaff_SF$scaffold)

plot(map_scaff_GL$X, map_scaff_GL$mean_cov, pch =20, cex = 0.8, xaxt="n",
     xlab="Scaffold", ylab="Average coverage", 
     main ="")#, ylim=c(0,100))

#abline(h = (quantile(map_scaff_GL$V3, probs = seq(0.995,1, 0.01), na.rm = T)), col= "grey40", lty=1)
points(map_scaff_LF$X, map_scaff_LF$mean_cov, pch = 2, cex = 0.8, col="red")
points(map_scaff_SF$X, map_scaff_SF$mean_cov, pch = 4, cex = 0.8, col="darkgreen")
points(map_scaff_LL$X, map_scaff_LL$mean_cov, pch = 18, cex = 0.8, col="blue")

legend("topleft", legend=c("German Lab", "London Field", "Sweden Field", "London Lab"), pch=c(20,2,4,18),
       col=c("black","red","darkgreen","blue"), cex=0.8, box.lty=1)
```

Summary for each population :

SF
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   0.00   10.22   17.08   23.84   19.25 5487.09 

LF
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.00   11.64   19.48   29.61   21.66 5963.07 

LL
Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.000    9.524   16.309   22.569   18.300 3986.625 

GL
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.00   15.04   25.21   37.73   28.02 8126.32 

<img src="order_length.png" width="600" height="500" style="background:none; border:none; box-shadow:none;">

As we can see, scaffold with extremely high mean coverage aren't the longer ones. Let's analyse it deeper:

```{r}
map_SF=read.table("map_scaff_SF.txt")
pvec <- seq(0.9,1,0.01)
quantile(map_SF$V3, pvec)
map_SF_highcov <- map_SF[map_SF$V3 > quantile(map_SF$V3, 0.99), ]

write.table(map_SF_highcov, file="/beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/SEP_MAP_UNMAP/map_scaff_SF_highcov.txt")

# Extraction of scaffold with high coverage (threshold 0.99) : 

map_SF=read.table("SF_cov_map.txt")
map_SF_highcov <- map_SF[ which(map_SF$V1=='NW_019392706.1' | map_SF$V1=='NW_019392715.1' | map_SF$V1=='NW_019392726.1' | map_SF$V1=='NW_019392782.1' | map_SF$V1=='NW_019392787.1' | map_SF$V1=='NW_019392930.1' | map_SF$V1=='NW_019393092.1' | map_SF$V1=='NW_019393097.1' | map_SF$V1=='NW_019393543.1' | map_SF$V1=='NW_019393765.1' | map_SF$V1=='NW_019393885.1' | map_SF$V1=='NW_019393980.1' | map_SF$V1=='NW_019394066.1' | map_SF$V1=='NW_019394087.1' | map_SF$V1=='NW_019394151.1' | map_SF$V1=='NC_030043.1'),]
write.table(map_SF_highcov, file="/beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/SEP_MAP_UNMAP/map_SF_highcov.txt")

# Threshold 0.995 :
map_SF_highcov <- map_SF[ which(map_SF$V1=='NW_019393765.1' | map_SF$V1=='NW_019393885.1' | map_SF$V1=='NW_019392930.1' | map_SF$V1=='NW_019392787.1' |  map_SF$V1=='NW_019392782.1' | map_SF$V1=='NW_019392715.1' | map_SF$V1=='NW_019392726.1' |    map_SF$V1=='NC_030043.1'),]
```

Extract scaffold sequences identified at the 0.995 threshold from reference genome and Blast with our own database, including several bacterian genome (Wolbachia_cimex_genome, g-proteobacteria_genome, clostridium_genome):

```
/beegfs/data/soft/samtools-1.9/bin/samtools faidx /beegfs/data/chaberkorn/PoolSeq_Clec/Ref_Clec/Cimex_lectularius.fna NW_019392715.1 > scaff_NW_019392715.fa

cat mitochondion_cimex.fasta wolbachia_cimex_genome.fna g-proteobacteria_genome.fna clostridium_genome.fna /beegfs/data/chaberkorn/PoolSeq_Clec/Ref_Clec/Cimex_lectularius.fna > data_cimex.fasta

/beegfs/data/chaberkorn/Tools/myconda/bin/makeblastdb -in data_cimex.fasta -dbtype nucl -out data_cimex 
# Output : 3 data_cimex.* file

# Blast also only with only bacterian genomes:

cat wolbachia_cimex_genome.fna g-proteobacteria_genome.fna clostridium_genome.fna > data_bacteria.fasta

/beegfs/data/chaberkorn/Tools/myconda/bin/makeblastdb -in data_bacteria.fasta -dbtype nucl -out data_bacteria 
/beegfs/data/chaberkorn/Tools/myconda/bin/blastn -query /beegfs/data/chaberkorn/PoolSeq_Clec/Genomes/ref_NC_030043.fa -db data_cimex -out scaff_030043_vs_data_cimex.blastn -outfmt 6 -max_target_seqs 5 -evalue 10e-1
```
Open "cimex_vs_data_TE.blastn" on Excel, convert in CSV:

```
data_TE <- read.csv("~/Desktop/CloudStation/THESE/WholeGenome PoolSeq/cimex_vs_data_TE.blastn.csv", sep=";",h=F)
colnames(data_TE) <- c("query_id","subject_id","identity","alignement_length","mismatches","gap_opens","qstart","qend","sstart","send","evalue","bitscore")

hist(log(data_TE$evalue)) # Artefact à 0
summary(data_TE$evalue)

# spliter par query pour mettre ensembles les scaffolds
liste=base::split(x=data_TE, f=data_TE$query_id) # 1300
length(liste)

# lapply
x=liste[[1]]
library(intervals)
i=Intervals(x[,c("qstart", "qend")])
union=interval_union(x = i)
as.data.frame(union)

liste_unions=lapply(X = liste, FUN = function(x){
  x=x[x$evalue<10^-10,]
  i=Intervals(x[,c("qstart", "qend")])
  union=interval_union(x = i)
  union=as.data.frame(union)
  return(union)
})

length(liste_unions)
liste_unions[1] # Premier élément = NC_030043.1 -> mitochondrial
liste_unions[2]

longueur=lapply(X = liste_unions, FUN=function(x){
  res=dim(x)[1]
  return(res)
})

head(longueur)

dataframe_bed=do.call(rbind.data.frame, liste_unions)
head(dataframe_bed)
noms_scaff=rep.int(names(liste_unions), times =longueur)
head(noms_scaff)
dataframe_bed$V3=noms_scaff
head(dataframe_bed)

dataframe_bed=dataframe_bed[, c(3,1,2)]
head(dataframe_bed)
dim(dataframe_bed) # 22914 lignes

hist(dataframe_bed$V2-dataframe_bed$V1) # Vérifier que tout est strictement >0 (end toujours > start)
dataframe_bed=dataframe_bed[-1,] # Enlever NC (mitochondrial)
colnames(dataframe_bed) <- c("query_id","qstart","qend")

write.table(dataframe_bed, file = "~/Desktop/CloudStation/THESE/WholeGenome PoolSeq/cimex_vs_data_TE.blastn.bed", row.names = F, quote=F)

# Ajouter la colonne subjecti_id
head(dataframe_bed)
dim(dataframe_bed) # 22913
head(data_TE)
dim(data_TE) # 78926 

library(dplyr)
test2 <- left_join(dataframe_bed, data_TE) %>% select(query_id,qstart,qend,subject_id)
head(test2)
dim(test2) # 26477 --> NA quand pas de TE correspondant à la région
View(data_TE)

# Marche pas : deux régions accolées en 1 unique --> ne retient aucune donnée de subject_id : traiter manuellement ?

# data_TE :
# NW_019392631.1	rnd-1_family-19#LINE/Penelope	4863	4920
#	NW_019392631.1	rnd-1_family-19#LINE/Penelope	4920	4993

# test2 :
# NW_019392631.1	4863	4993	NA

# Contre-exemple : deux régions superposées --> ne retient que la première ligne

# data_TE :
# NW_019392631.1	rnd-1_family-58#LINE/LOA	10269	11010	
#	NW_019392631.1	rnd-6_family-1023#LINE/LOA	10277	11010

# test2 :
# NW_019392631.1	10269	11010	rnd-1_family-58#LINE/LOA
```

## Genetic differenciation of populations

The goal is to understand what differenciates the four strains of *Cimex lectularius* PoolSeq samples: London Lab, London Field, German Lab, Sweden Field. 
Which are the most similar? What seems to differentiate them?

hypopthèse : certaines R/ d'autres S (donner infos site cimesxtore)


We will have to download a few soft.

### Install tools

PoPoolation & PoPoolation2 :
``` 
chaberkorn@pbil-deb:/beegfs/data/chaberkorn/Tools/popoolation_1.2.2
chaberkorn@pbil-deb:/beegfs/data/chaberkorn/Tools/popoolation2_1201
```

Poolfstat :
``` 
install.packages("poolfstat")
```


### Detecting SNPs

A SNP (Single Nucleotide Polymorphism) is a ponctual mutation that may be associated with candidate regions for resistance.
The output file 'all.sync' from PoPoolation can be used with Poolfstat to compute the SNPs:  

``` 
library(poolfstat)

# wc -l all.sync = 509,820,671
# Estimate time - 0.36s/Mi lines processed -> 510*0.36 ~ 184 min soit 3h
# Default parameters - start at 13h45, end at ~17h (189.86 min)

pooldata = popsync2pooldata(sync.file = "all.sync", poolsizes = rep(30,4), 
                            poolnames = c("GL","LF","SF","LL"), min.rc = 1, min.cov.per.pool = -1,
                            max.cov.per.pool = 1e+06, min.maf = 0.01, noindel = TRUE,
                            nlines.per.readblock = 1e+06, nthreads = 1)
# 509.8207 millions lines processed in 189.86  min.;  10 139 943 SNPs found for 4 Pools

# Do subset of pooldata :
pooldata_sub <- pooldata.subset(pooldata, pool.index = c(1,2,3,4), min.cov.per.pool = 10, max.cov.per.pool = 200, min.maf = 0.05)
# pool.index = indexes of the pools (at least two), that should be selected to create the new pooldata objec
# min.maf correspond to Minimal allowed Minor Allele Frequency -> min-count of 5 w/ PoPoolation2
# With same parameters, 5 722 762 SNPs for 4 Pools, compared to 5 837 216 w/ PoPoolation2)
``` 

Perform PCA:

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

# Different (?) from : SNP_biall$fq_SF <- ifelse(SNP_biall$ref==SNP_biall$SF1, SNP_biall$SF_ref/SNP_biall$SF_tot, SNP_biall$SF_alt/SNP_biall$SF_tot) 
  
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

<img src="ACP_poolfstat.png" width="600" height="500" style="background:none; border:none; box-shadow:none;">

PCA allows us to detect the SNPs distinguishing the 4 strains: London Field and London Lab seems poorly genetically differenciated. The genetic differenciation appears to be related to the geographical repartition. 



### Compute FST

``` 
PairwiseFST = computePairwiseFSTmatrix(pooldata_sub, method = "Anova", min.cov.per.pool = -1, 
                                     max.cov.per.pool = 1e+06, min.maf = -1, output.snp.values = FALSE)
# output.snp.values : If TRUE, provide SNP-specific pairwise FST for each comparisons (may lead to a huge result object if the number of pools and/or SNPs is large)

PairwiseFST_all = na.omit(computePairwiseFSTmatrix(pooldata_sub, method = "Anova", min.cov.per.pool = -1, 
                                     max.cov.per.pool = 1e+06, min.maf = -1, output.snp.values = TRUE))
``` 

Visualize FST:

``` 
PairwiseFST_all$PairwiseSnpFST
# Ordre : GL_vs_LF, GL_vs_SF, GL_vs_LL, LF_vs_SF, LF_vs_LL, SF_vs_LL

mean(PairwiseFST_all$PairwiseSnpFST[,5], na.rm=T)
summary(PairwiseFST_all$PairwiseSnpFST)

boxplot(PairwiseFST_all$PairwiseSnpFST[,1], PairwiseFST_all$PairwiseSnpFST[,2], PairwiseFST_all$PairwiseSnpFST[,3], PairwiseFST_all$PairwiseSnpFST[,4], PairwiseFST_all$PairwiseSnpFST[,5], PairwiseFST_all$PairwiseSnpFST[,6], 
        border=c("#F8766D","#D89000","#72B000","#00C19C","#00B0F6","#CF78FF"), 
        xlab="Strains comparisons", 
        ylab="FSTs distribution",
        names=c("GL_vs_LF","GL_vs_SF","GL_vs_LL","LF_vs_SF","LF_vs_LL","SF_vs_LL"))
``` 

<img src="boxplot_poolfstat.png" width="600" height="400" style="background:none; border:none; box-shadow:none;">

Identify outliers FST:

``` 
FST_tab_LG$GL_vs_LF <- as.numeric(FST_tab_LG$GL_vs_LF)
FST_tab_LG$GL_vs_SF <- as.numeric(FST_tab_LG$GL_vs_SF)
FST_tab_LG$LF_vs_LL <- as.numeric(FST_tab_LG$LF_vs_LL)
FST_tab_LG$SF_vs_LL <- as.numeric(FST_tab_LG$SF_vs_LL)
        
FST_tab_LG$Colour="black"
FST_tab_LG$Colour[(FST_tab_LG$GL_vs_LF > quantile(FST_tab_LG$GL_vs_LF, 0.99, na.rm=T) 
                   & FST_tab_LG$GL_vs_SF > quantile(FST_tab_LG$GL_vs_SF, 0.99, na.rm=T)
                   & FST_tab_LG$LF_vs_LL > quantile(FST_tab_LG$LF_vs_LL, 0.99, na.rm=T)
                   & FST_tab_LG$SF_vs_LL > quantile(FST_tab_LG$SF_vs_LL, 0.99, na.rm=T))]="red"

FST_tab_LG$Colour[(FST_tab_LG$GL_vs_LF > quantile(FST_tab_LG$GL_vs_LF, 0.99, na.rm=T) 
                   & FST_tab_LG$GL_vs_SF > quantile(FST_tab_LG$GL_vs_SF, 0.99, na.rm=T)
                   & FST_tab_LG$LF_vs_LL > quantile(FST_tab_LG$LF_vs_LL, 0.99, na.rm=T)
                   & FST_tab_LG$SF_vs_LL > quantile(FST_tab_LG$SF_vs_LL, 0.99, na.rm=T))]="red"

length((FST_tab_LG$GL_vs_LF > quantile(FST_tab_LG$GL_vs_LF, 0.99, na.rm=T) & FST_tab_LG$GL_vs_SF > quantile(FST_tab_LG$GL_vs_SF, 0.99, na.rm=T))[(FST_tab_LG$GL_vs_LF > quantile(FST_tab_LG$GL_vs_LF, 0.99, na.rm=T) & FST_tab_LG$GL_vs_SF > quantile(FST_tab_LG$GL_vs_SF, 0.99, na.rm=T))==TRUE])

#  Missing values and NaN are not allowed if 'na.rm' is FALSE (default value)

library(dplyr)
FST_tab_LG <- FST_tab_LG %>% arrange(Colour)

# Define subset for each chromosome :

chr1 <- subset(FST_tab_LG, LG == 1)
chr2 <- subset(FST_tab_LG, LG == 2)
chr3 <- subset(FST_tab_LG, LG == 3)
chr4 <- subset(FST_tab_LG, LG == 4)
chr5 <- subset(FST_tab_LG, LG == 5)
chr6 <- subset(FST_tab_LG, LG == 6)
chr7 <- subset(FST_tab_LG, LG == 7)
chr8 <- subset(FST_tab_LG, LG == 8)
chr9 <- subset(FST_tab_LG, LG == 9)
chr10 <- subset(FST_tab_LG, LG == 10)
chr11 <- subset(FST_tab_LG, LG == 11)
chr12 <- subset(FST_tab_LG, LG == 12)
chr13 <- subset(FST_tab_LG, LG == 13)
chr14 <- subset(FST_tab_LG, LG == 14)

# Open a script and copy paste for each combination :

par(mfrow=c(3,5))
plot(x=chr1$position,y=chr1$LF_vs_LL, ylim = c(0,1), ylab="", xlab="", pch=20, main="Chr 1", col=chr1$Colour)
plot(x=chr2$position,y=chr2$LF_vs_LL, ylim = c(0,1), ylab="", xlab="", pch=20, main="Chr 2", col=chr2$Colour)
plot(x=chr3$position,y=chr3$LF_vs_LL, ylim = c(0,1), ylab="", xlab="", pch=20, main="Chr 3", col=chr3$Colour)
plot(x=chr4$position,y=chr4$LF_vs_LL, ylim = c(0,1), ylab="", xlab="", pch=20, main="Chr 4", col=chr4$Colour)
plot(x=chr5$position,y=chr5$LF_vs_LL, ylim = c(0,1), ylab="", xlab="", pch=20, main="Chr 5", col=chr5$Colour)
plot(x=chr6$position,y=chr6$LF_vs_LL, ylim = c(0,1), ylab="", xlab="", pch=20, main="Chr 6", col=chr6$Colour)
plot(x=chr7$position,y=chr7$LF_vs_LL, ylim = c(0,1), ylab="", xlab="", pch=20, main="Chr 7", col=chr7$Colour)
plot(x=chr8$position,y=chr8$LF_vs_LL, ylim = c(0,1), ylab="", xlab="", pch=20, main="Chr 8", col=chr8$Colour)
plot(x=chr9$position,y=chr9$LF_vs_LL, ylim = c(0,1), ylab="", xlab="", pch=20, main="Chr 9", col=chr9$Colour)
plot(x=chr10$position,y=chr10$LF_vs_LL, ylim = c(0,1), ylab="", xlab="", pch=20, main="Chr 10", col=chr10$Colour)
plot(x=chr11$position,y=chr11$LF_vs_LL, ylim = c(0,1), ylab="", xlab="", pch=20, main="Chr 11", col=chr11$Colour)
plot(x=chr12$position,y=chr12$LF_vs_LL, ylim = c(0,1), ylab="", xlab="", pch=20, main="Chr 12", col=chr12$Colour)
plot(x=chr13$position,y=chr13$LF_vs_LL, ylim = c(0,1), ylab="", xlab="", pch=20, main="Chr 13", col=chr13$Colour)
plot(x=chr14$position,y=chr14$LF_vs_LL, ylim = c(0,1), ylab="", xlab="", pch=20, main="Chr 14", col=chr14$Colour)
``` 
![Image 5](FST_LF_LL_poolfstat.png)

``` 
poolfstat_high_FST_LG <- subset(FST_tab_LG, Colour=="red")
write.table(poolfstat_high_FST_LG, file="poolfstat_high_FST_LG.txt", sep=",")

# How many high FST in each chr ?

chr1_high <- subset(FST_tab_LG, LG == 1 & Colour == "red")
chr2_high <- subset(FST_tab_LG, LG == 2 & Colour == "red")
chr3_high <- subset(FST_tab_LG, LG == 3 & Colour == "red")
chr4_high <- subset(FST_tab_LG, LG == 4 & Colour == "red")
chr5_high <- subset(FST_tab_LG, LG == 5 & Colour == "red")
chr6_high <- subset(FST_tab_LG, LG == 6 & Colour == "red")
chr7_high <- subset(FST_tab_LG, LG == 7 & Colour == "red")
chr8_high <- subset(FST_tab_LG, LG == 8 & Colour == "red")
chr9_high <- subset(FST_tab_LG, LG == 9 & Colour == "red")
chr10_high <- subset(FST_tab_LG, LG == 10 & Colour == "red")
chr11_high <- subset(FST_tab_LG, LG == 11 & Colour == "red")
chr12_high <- subset(FST_tab_LG, LG == 12 & Colour == "red")
chr13_high <- subset(FST_tab_LG, LG == 13 & Colour == "red")
chr14_high <- subset(FST_tab_LG, LG == 14 & Colour == "red")

# Find associate genes

scaff_genes=read.table("scaff_genes.txt") 
scaff_genes$seqid <- substring(scaff_genes$seqid,1,12)

# A faire sur le cluster :
# scp /Users/chloe/Documents/Cluster/poolfstat_high_FST_LG.txt chaberkorn@pbil-gates.univ-lyon1.fr:/beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/SEP_MAP_UNMAP/
# cd /beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/SEP_MAP_UNMAP/
# /beegfs/data/soft/R-3.5.2/bin/R

scaff_genes=read.table("scaff_genes.txt") 

library(dplyr)
library(tidyr)

newtest_genes <- 
  scaff_genes %>% 
  rowwise %>% 
  mutate(pos = paste(seq(start, 
                         end), 
                     collapse = ",")) %>%
  tidyr::separate_rows(pos) %>%
  mutate(pos = as.integer(pos))

newtest_genes<-subset(newtest_genes, select=-c(source,type,score,strand,phase,ID,
                                               Dbxref,Name,gbkey,gene_biotype))
colnames(newtest_genes) <- c("scaffold", "start","end","gene","length","position")
newtest_genes$scaffold <- substring(newtest_genes$scaffold,1,12)

high_FST=read.table("poolfstat_high_FST_LG.txt", header=T, sep=",") 
poolfstat_genes_high_FST_LG <- right_join(newtest_genes, high_FST) # Joining, by = c("scaffold", "position")

write.table(newtest_genes, file="/beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/SEP_MAP_UNMAP/newtest_genes.txt", sep=",")
write.table(poolfstat_genes_high_FST_LG, file="/beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/SEP_MAP_UNMAP/poolfstat_genes_high_FST_LG.txt", sep=",")

# Back to R computer session

poolfstat_genes_high_FST_LG <- read.table(file="~/Desktop/CloudStation/THESE/WholeGenome PoolSeq/SNPs/poolfstat_genes_high_FST_LG.txt", sep=",", header=T)
View(poolfstat_genes_high_FST_LG)

# IF we have gene ID, merge with scaff_cds table to have product column

scaff_cds=read.table("~/Desktop/CloudStation/THESE/CNV/scaff_cds.txt") 
head(scaff_cds)

library(tidyverse)
scaff_cds <- scaff_cds %>%
  select(gene, product)
View(scaff_cds)

poolfstat_cds_high_FST_LG <- left_join(poolfstat_genes_high_FST_LG,scaff_cds, by=c("gene"))
View(high_FST_cds)

uniq_poolfstat_cds_high_FST_LG <- poolfstat_cds_high_FST_LG[!duplicated(poolfstat_cds_high_FST_LG[,c('start','position')]),]
uniq_poolfstat_cds_high_FST_LG$X <- 1:length(uniq_poolfstat_cds_high_FST_LG$X)
View(uniq_poolfstat_cds_high_FST_LG)
```

UGPMA can be compute on FST data :

``` 
library("phangorn")
matrix_upgma_auto <- PairwiseFST_all$PairwiseFSTmatrix
# Or manually:
matrix_upgma=matrix(nrow=4,ncol=4) # pour LG, sinon changer 
colnames(matrix_upgma)=c("GL","LF","SF","LL")
rownames(matrix_upgma)=c("GL","LF","SF","LL")
matrix_upgma[c(1,2,3,4),c(1,2,3,4)]=0

# Order : GL_vs_LF, GL_vs_SF, GL_vs_LL, LF_vs_SF, LF_vs_LL, SF_vs_LL
matrix_upgma[2,1]=mean(PairwiseFST_all$PairwiseSnpFST[,1], na.rm=T)
matrix_upgma[1,2]=mean(PairwiseFST_all$PairwiseSnpFST[,1], na.rm=T)

matrix_upgma[3,1]=mean(PairwiseFST_all$PairwiseSnpFST[,2], na.rm=T)
matrix_upgma[1,3]=mean(PairwiseFST_all$PairwiseSnpFST[,2], na.rm=T)

matrix_upgma[4,1]=mean(PairwiseFST_all$PairwiseSnpFST[,3], na.rm=T)
matrix_upgma[1,4]=mean(PairwiseFST_all$PairwiseSnpFST[,3], na.rm=T)

matrix_upgma[2,3]=mean(PairwiseFST_all$PairwiseSnpFST[,4], na.rm=T)
matrix_upgma[3,2]=mean(PairwiseFST_all$PairwiseSnpFST[,4], na.rm=T)

matrix_upgma[2,4]=mean(PairwiseFST_all$PairwiseSnpFST[,5], na.rm=T)
matrix_upgma[4,2]=mean(PairwiseFST_all$PairwiseSnpFST[,5], na.rm=T)

matrix_upgma[3,4]=mean(PairwiseFST_all$PairwiseSnpFST[,6], na.rm=T)
matrix_upgma[4,3]=mean(PairwiseFST_all$PairwiseSnpFST[,6], na.rm=T)

data_upgma_auto <- upgma(matrix_upgma_auto, method = "average")
data_upgma <- upgma(matrix_upgma, method = "average")

plot(data_upgma)
plot(data_upgma_auto)
``` 
![Image 7](.png)


### Estimate genetic polymorphism



Convert to use on BayPass :

```{r}
setwd("~/Desktop/CloudStation/THESE/WholeGenome PoolSeq/SNPs/Poolfstat")

pooldata2genobaypass(pooldata,
                     writing.dir = getwd(), # directory where to create the files
                     prefix = "poolfstatdata_110221", # prefix used for output file names
                     subsamplesize = -1) # all SNPs are considered in the output
                     # subsamplingmethod = "thinning")

pooldata2genobaypass(pooldata_sub,
                     writing.dir = getwd(), # directory where to create the files
                     prefix = "poolfstatdata_sub_110221", # prefix used for output file names
                     subsamplesize = -1) # all SNPs are considered in the output
                     # subsamplingmethod = "thinning")
```









