#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
#args[1]<-"/Users/mcallister/Desktop/test_out/blast_results" #working directory
#args[2]<-326 #expected amplicon size (excluding primers)
########################################
library(dplyr)
setwd(as.character(args[1]))

all_results <- read.delim("ASV_blastn_nt.btab", header=FALSE, stringsAsFactors=FALSE)
colnames(all_results) <- c("ASV", "percent", "length", "taxid")
all_results_trimLength <- all_results[all_results$length>=(0.9*as.numeric(args[2])), ]
all_results_trimLength$correction <- all_results_trimLength$percent/100 *all_results_trimLength$length
all_results_trimLength_filter <- all_results_trimLength %>% group_by(ASV) %>% filter(percent == max(percent))
all_results_bestHit <- all_results_trimLength_filter %>% group_by(ASV)%>%summarise_all(funs(toString(unique(.))))
write.table(all_results_bestHit, file='ASV_blastn_nt_formatted.txt', sep='\t', quote=FALSE, row.names=FALSE)