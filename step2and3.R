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