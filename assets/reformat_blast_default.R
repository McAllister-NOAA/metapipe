#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
args[1]<-"/Users/mcallister/Desktop/test_out/blast_results" #working directory
args[2]<-326 #expected amplicon size (excluding primers)

########################################
library(dplyr)
setwd(as.character(args[1]))

all_results <- read.delim("ASV_blastn_nt.btab", header=FALSE, stringsAsFactors=FALSE)
ignoreMostEnvSeqs <- read.delim("ASV_blastn_nt_ignoreMostEnvSeqs_leavesinUnclassified.btab", header=FALSE, stringsAsFactors=FALSE)
ignoreAllEnvSeqs <- read.delim("ASV_blastn_nt_ignoreAllEnvSeqs_removesUnclassified.btab", header=FALSE, stringsAsFactors=FALSE)

colnames(all_results) <- c("ASV", "percent", "length", "taxid")
colnames(ignoreMostEnvSeqs) <- c("ASV", "percent", "length", "taxid")
colnames(ignoreAllEnvSeqs) <- c("ASV", "percent", "length", "taxid")


all_results_trimLength <- all_results[all_results$length>=(0.9*as.numeric(args[2])), ]
ignoreMostEnvSeqs_trimLength <- ignoreMostEnvSeqs[ignoreMostEnvSeqs$length>=(0.9*as.numeric(args[2])), ]
ignoreAllEnvSeqs_trimLength <- ignoreAllEnvSeqs[ignoreAllEnvSeqs$length>=(0.9*as.numeric(args[2])), ]

all_results_trimLength$correction <- all_results_trimLength$percent/100 *all_results_trimLength$length
ignoreMostEnvSeqs_trimLength$correction <- ignoreMostEnvSeqs_trimLength$percent/100 *ignoreMostEnvSeqs_trimLength$length
ignoreAllEnvSeqs_trimLength$correction <- ignoreAllEnvSeqs_trimLength$percent/100 *ignoreAllEnvSeqs_trimLength$length

all_results_trimLength_filter <- all_results_trimLength %>% group_by(ASV) %>% filter(percent == max(percent))
ignoreMostEnvSeqs_trimLength_filter <- ignoreMostEnvSeqs_trimLength %>% group_by(ASV) %>% filter(percent == max(percent))
ignoreAllEnvSeqs_trimLength_filter <- ignoreAllEnvSeqs_trimLength %>% group_by(ASV) %>% filter(percent == max(percent))


all_results_bestHit <- all_results_trimLength_filter %>% group_by(ASV)%>%summarise_all(funs(toString(unique(.))))
ignoreMostEnvSeqs_bestHit <- ignoreMostEnvSeqs_trimLength_filter %>% group_by(ASV)%>%summarise_all(funs(toString(unique(.))))
ignoreAllEnvSeqs_bestHit <- ignoreAllEnvSeqs_trimLength_filter %>% group_by(ASV)%>%summarise_all(funs(toString(unique(.))))

write.table(all_results_bestHit, file='ASV_blastn_nt_formatted.txt', sep='\t', quote=FALSE, row.names=FALSE)
write.table(ignoreMostEnvSeqs_bestHit, file='ASV_blastn_nt_ignoreMostEnvSeqs_leavesinUnclassified_formatted.txt', sep='\t', quote=FALSE, row.names=FALSE)
write.table(ignoreAllEnvSeqs_bestHit, file='ASV_blastn_nt_ignoreAllEnvSeqs_removesUnclassified_formatted.txt', sep='\t', quote=FALSE, row.names=FALSE)
