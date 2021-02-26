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

Represent coverages on R:

```
setwd("~/Desktop/CloudStation/THESE/WholeGenome PoolSeq/Analyse coverage")

map_scaff_GL=read.table("map_scaff_GL.txt")
map_scaff_LL=read.table("map_scaff_LL.txt")
map_scaff_LF=read.table("map_scaff_LF.txt")
map_scaff_SF=read.table("map_scaff_SF.txt")

library(ggplot2)

head(map_scaff_GL)

GL_highcov <- map_scaff_GL[ which(map_scaff_GL$scaffold=='NW_019392715.1' | map_scaff_GL$scaffold=='NW_019392726.1' | map_scaff_GL$scaffold=='NW_019392782.1' | map_scaff_GL$scaffold=='NW_019392787.1' | map_scaff_GL$scaffold=='NW_019392930.1' |  map_scaff_GL$scaffold=='NW_019393765.1' | map_scaff_GL$scaffold=='NW_019393885.1' | map_scaff_GL$scaffold=='NC_030043.1'),]

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

![Image 2](bedbugs.png)
![Image 3](bedbugs.png)





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

Analyser les hautes couvertures :

```{r}
map_SF=read.table("map_scaff_SF.txt")

pvec <- seq(0,1,0.1)
quantile(map_SF$V3, pvec)
pvec <- seq(0.9,1,0.01)
quantile(map_SF$V3, pvec)
map_SF_highcov <- map_SF[map_SF$V3 > quantile(map_SF$V3, 0.99), ]

write.table(map_SF_highcov, file="/beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/SEP_MAP_UNMAP/map_scaff_SF_highcov.txt")

# On veut extraire les scaffold ayant une haute couverture du tableau global (seuil 0.99):

map_SF=read.table("SF_cov_map.txt")

map_SF_highcov <- map_SF[ which(map_SF$V1=='NW_019392706.1' | map_SF$V1=='NW_019392715.1' | map_SF$V1=='NW_019392726.1' | map_SF$V1=='NW_019392782.1' | map_SF$V1=='NW_019392787.1' | map_SF$V1=='NW_019392930.1' | map_SF$V1=='NW_019393092.1' | map_SF$V1=='NW_019393097.1' | map_SF$V1=='NW_019393543.1' | map_SF$V1=='NW_019393765.1' | map_SF$V1=='NW_019393885.1' | map_SF$V1=='NW_019393980.1' | map_SF$V1=='NW_019394066.1' | map_SF$V1=='NW_019394087.1' | map_SF$V1=='NW_019394151.1' | map_SF$V1=='NC_030043.1'),]
write.table(map_SF_highcov, file="/beegfs/data/chaberkorn/PoolSeq_Clec/Mapped/SEP_MAP_UNMAP/map_SF_highcov.txt")

# Seuil 0.995 :
map_SF_highcov <- map_SF[ which(map_SF$V1=='NW_019393765.1' | map_SF$V1=='NW_019393885.1' | map_SF$V1=='NW_019392930.1' | map_SF$V1=='NW_019392787.1' |  map_SF$V1=='NW_019392782.1' | map_SF$V1=='NW_019392715.1' | map_SF$V1=='NW_019392726.1' |    map_SF$V1=='NC_030043.1'),]

scaff_1 <- map_SF_highcov[ which(map_SF_highcov$V1=='NW_019392715.1'),]
```

Regarder manuellement dans les scaffolds identifies au seuil 0.995 ou sont les couvertures >1000, et exporter ensuite les sequences de ces positions en les extrayant du genome de reference :

```{bash}
/beegfs/data/soft/samtools-1.9/bin/samtools faidx /beegfs/data/chaberkorn/PoolSeq_Clec/Ref_Clec/Cimex_lectularius.fna NW_019392715.1:376100-376700 > ref_NW_019392715.fai

/beegfs/data/soft/samtools-1.9/bin/samtools faidx /beegfs/data/chaberkorn/PoolSeq_Clec/Ref_Clec/Cimex_lectularius.fna NW_019392726.1:966100-980300 > ref_NW_019392726.fai

/beegfs/data/soft/samtools-1.9/bin/samtools faidx /beegfs/data/chaberkorn/PoolSeq_Clec/Ref_Clec/Cimex_lectularius.fna NW_019392782.1:1283600-1337700 > ref_NW_019392782.fai

/beegfs/data/soft/samtools-1.9/bin/samtools faidx /beegfs/data/chaberkorn/PoolSeq_Clec/Ref_Clec/Cimex_lectularius.fna NW_019392787.1:450700-450900 > ref_NW_019392787.fai

/beegfs/data/soft/samtools-1.9/bin/samtools faidx /beegfs/data/chaberkorn/PoolSeq_Clec/Ref_Clec/Cimex_lectularius.fna NW_019392930.1:195300-195600  > ref_NW_019392930.fai

/beegfs/data/soft/samtools-1.9/bin/samtools faidx /beegfs/data/chaberkorn/PoolSeq_Clec/Ref_Clec/Cimex_lectularius.fna NW_019393765.1:480-750 > ref_NW_019393765.fai

/beegfs/data/soft/samtools-1.9/bin/samtools faidx /beegfs/data/chaberkorn/PoolSeq_Clec/Ref_Clec/Cimex_lectularius.fna NW_019393885.1:160- > ref_NW_019393885.fai

# Ou pour blaster :
/beegfs/data/soft/samtools-1.9/bin/samtools faidx /beegfs/data/chaberkorn/PoolSeq_Clec/Ref_Clec/Cimex_lectularius.fna NC_030043.1 > scaff_NC_030043.fa
```

Utiliser Blast en ligne de commande en creeant notre propre base de donnes pour identifier la nature de ces scaffold hautement couverts :

```{bash}
cat mitochondion_cimex.fasta wolbachia_cimex_genome.fna g-proteobacteria_genome.fna clostridium_genome.fna /beegfs/data/chaberkorn/PoolSeq_Clec/Ref_Clec/Cimex_lectularius.fna > data_cimex.fasta

/beegfs/data/chaberkorn/Tools/myconda/bin/makeblastdb -in data_cimex.fasta -dbtype nucl -out data_cimex 
# Output : 3 fichiers data_cimex.*

# Ou avec uniquement genome bacterien :

cat wolbachia_cimex_genome.fna g-proteobacteria_genome.fna clostridium_genome.fna > data_bacteria.fasta

/beegfs/data/chaberkorn/Tools/myconda/bin/makeblastdb -in data_bacteria.fasta -dbtype nucl -out data_bacteria 

# Puis blaster :

/beegfs/data/chaberkorn/Tools/myconda/bin/blastn -query /beegfs/data/chaberkorn/PoolSeq_Clec/Genomes/ref_NC_030043.fa -db data_cimex -out scaff_030043_vs_data_cimex.blastn -outfmt 6 -max_target_seqs 5 -evalue 10e-1
```

Ouvrir "cimex_vs_data_TE.blastn" sur Excel, convertir en CSV

```{r}
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
















Poolation & PoPoolation2 :

``` 
chaberkorn@pbil-deb:/beegfs/data/chaberkorn/Tools/popoolation_1.2.2
chaberkorn@pbil-deb:/beegfs/data/chaberkorn/Tools/popoolation2_1201
```
