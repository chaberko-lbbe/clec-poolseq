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
  # file1="/your-path/PoolSeq_Clec/CNV/result-cov-random/random-windows-10kb_covLL_ok.txt" 
  # file2="/your-path/PoolSeq_Clec/CNV/result-cov-random/random-windows-10kb_covLF_ok.txt" 
  # output_file_1="/your-path/PoolSeq_Clec/CNV/step4pt5/data/RDR_95/quantiles_10kb.txt"
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

LL=read.table(file1, sep="\t", header=F)
LF=read.table(file2, sep="\t", header=F)

##################################################################################################################
# Compute ratio LL/LF and its quantiles
##################################################################################################################

colnames(LL) <- c("rname","startpos","endpos","sumsize","sumcov")
colnames(LF) <- c("rname","startpos","endpos","sumsize","sumcov")

LL$ratio <- (LL$sumcov/LL$sumsize*100)/(LF$sumcov/LF$sumsize*100)
q05 <- quantile(LL$ratio,prob=0.05,na.rm=T)[[1]]
q95 <- quantile(LL$ratio,prob=0.95,na.rm=T)[[1]]
q10 <- quantile(LL$ratio,prob=0.10,na.rm=T)[[1]]
q90 <- quantile(LL$ratio,prob=0.90,na.rm=T)[[1]]

quantiles <- c(q05,q10,q90,q95)

##################################################################################################################
# Write to disk
##################################################################################################################

write.table(x = quantiles, file = output_file_1, quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")

