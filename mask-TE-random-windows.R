#!/usr/bin/env Rscript

##################################################################################################################
# parse arguments given to the script
##################################################################################################################

args = commandArgs(trailingOnly=TRUE)

# test if there is the required number of argument: if not, return an error
if (length(args)!=3) {
  stop("Three arguments must be supplied", call.=FALSE)
}
if (length(args)==3) {
  # file1="/your-path/PoolSeq_Clec/CNV/random-windows-500.txt.bed" 
  # file2="/your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/random-windows-500.txt_covLL.txt" 
  # output_file_1="/your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/random-windows-500_covLL_ok.txt"
  file1=args[1]
  file2=args[2]
  output_file_1=args[3] 
}

##################################################################################################################
#  import data 
##################################################################################################################

# Load librairies
suppressPackageStartupMessages({
  library(dplyr)
})

original_bed=read.table(file1, sep="\t", header=F)
cov_bed=read.table(file2, sep="\t", header=F)

##################################################################################################################
# Merge the two bed files back together:
##################################################################################################################

colnames(original_bed) <- c("scaffold","start","end") # 1,000
colnames(cov_bed) <- c("scaffold","start","end","cov") # 1,095
# cov =  total read base count (i.e. the sum of per base read depths) 

cov_bed$size <- cov_bed$end-cov_bed$start

true_bed <- merge(cov_bed, original_bed, by=c("scaffold"), all.x=T, all.y=T)
dim(true_bed) # 5029

true_cov_bed <- true_bed %>% 
  group_by(scaffold,start.y,end.y) %>% 
  summarize(sumsize = sum(size[start.x >= start.y & end.x <= end.y]),
            sumcov = sum(cov[start.x >= start.y & end.x <= end.y]))

##################################################################################################################
# write to disk
##################################################################################################################

write.table(x = true_cov_bed, file = output_file_1, quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")

