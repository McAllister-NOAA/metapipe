#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
args[1]<-"/Users/mcallister/Desktop/test_figs" #FIGURE OUT directory
args[2]<-"/Users/mcallister/Desktop/chris_test/CP_all_out/ASV2Taxonomy/ASVs_counts_NOUNKNOWNS.tsv" #ASV count table
args[3]<-"/Users/mcallister/Desktop/chris_test/CP_all_out/ASV2Taxonomy/CP_all_out_asvTaxonomyTable_NOUNKNOWNS.txt" #ASV taxonomy table
args[4]<-"/Users/mcallister/Desktop/chris_test/CP_all_out/sample_metadata_forR.txt" #sample metadata file
args[5]<-5 #percent filter
args[6]<-TRUE #replicateFlag
args[7]<-FALSE #sitelabelFlag
########################################
library("ggplot2")
library("dplyr")
library("phyloseq")
library("ggpubr")
library("ggalt")
#library("plyr")

setwd(as.character(args[1]))
theme_set(theme_bw())

asv_count <- read.delim(as.character(args[2]), header=TRUE, stringsAsFactors=FALSE)
row.names(asv_count) <- asv_count$x
asv_count <- asv_count %>% select(-x)
asv_count_mat <- as.matrix(asv_count)

asv_taxonomy <- read.delim(as.character(args[3]), header=TRUE, stringsAsFactors=FALSE)
row.names(asv_taxonomy) <- asv_taxonomy$ASV
asv_taxonomy <- asv_taxonomy %>% select(-ASV)
asv_taxonomy_mat <- as.matrix(asv_taxonomy)

sample_metadata <- read.delim(as.character(args[4]), header=TRUE, stringsAsFactors=TRUE)
row.names(sample_metadata) <- sample_metadata$Sample
sample_metadata <- sample_metadata %>% select(-Sample)

ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
TAX = tax_table(asv_taxonomy_mat)
samples = sample_data(sample_metadata)

phylo_combined <- phyloseq(ASV, TAX, samples)

total = median(sample_sums(phylo_combined))
standf = function(x, t=total) round(t * (x / sum(x)))
phylo_combined_norm <- transform_sample_counts(phylo_combined, standf)

##BARGRAPHS (NO COLLAPSE) - NON-NORMALIZED
setwd(paste0(as.character(args[1]), "/nonnormalized", sep = ""))

for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  bar_plot <- plot_bar(phylo_combined, fill = i) +
    geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack")
  bar_plot_legend <- get_legend(bar_plot)
  bar_plot <- bar_plot + theme(legend.position='none')
  pdf(file=paste0('barplot_nonnormalized_',i,'_legend.pdf', sep = ""), width = 22, height = 17)
  print(as_ggplot(bar_plot_legend))
  dev.off()
  pdf(file=paste0('barplot_nonnormalized_',i,'.pdf', sep = ""), width = 11, height = 8.5)
  print(bar_plot)
  dev.off()
}

##BARGRAPHS (NO COLLAPSE) - NORMALIZED
setwd(paste0(as.character(args[1]), "/normalized", sep = ""))

for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  bar_plot <- plot_bar(phylo_combined_norm, fill = i) +
    geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack")
  bar_plot_legend <- get_legend(bar_plot)
  bar_plot <- bar_plot + theme(legend.position='none')
  pdf(file=paste0('barplot_normalized_',i,'_legend.pdf', sep = ""), width = 22, height = 17)
  print(as_ggplot(bar_plot_legend))
  dev.off()
  pdf(file=paste0('barplot_normalized_',i,'.pdf', sep = ""), width = 11, height = 8.5)
  print(bar_plot)
  dev.off()
}

#Filtering
filter_percent <- as.numeric(args[5])
# 
# asv_count <- colSums(asv_count[-1, ])
# 
# abundance_dat$readSum <- rowSums(abundance_dat[, -1])
# abundance_dat <- abundance_dat %>% mutate_at(vars(-Sample, -readSum), funs(100*./readSum))
# abundance_dat <- abundance_dat %>% select(-readSum)


##GROUP BY Replicates or Sites
replicateFlag <- as.logical(args[6])
sitelabelFlag <- as.logical(args[7])


if (sitelabelFlag == TRUE) { #Will expect "sites" column
  setwd(paste0(as.character(args[1]), "/nonnormalized", sep = ""))
  phylo_combined_reps <- merge_samples(phylo_combined, "sites")
  for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    bar_plot <- plot_bar(phylo_combined_reps, fill = i) +
      geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack")
    bar_plot_legend <- get_legend(bar_plot)
    bar_plot <- bar_plot + theme(legend.position='none')
    pdf(file=paste0('barplot_site_groups_nonnormalized_',i,'_legend.pdf', sep = ""), width = 22, height = 17)
    print(as_ggplot(bar_plot_legend))
    dev.off()
    pdf(file=paste0('barplot_site_groups_nonnormalized_',i,'.pdf', sep = ""), width = 11, height = 8.5)
    print(bar_plot)
    dev.off()
  }
  setwd(paste0(as.character(args[1]), "/normalized", sep = ""))
  newtotal = median(sample_sums(phylo_combined_reps))
  standf = function(x, t=total) round(t * (x / sum(x)))
  phylo_combined_norm_reps <- transform_sample_counts(phylo_combined_reps, standf)
  for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    bar_plot <- plot_bar(phylo_combined_norm_reps, fill = i) +
      geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack")
    bar_plot_legend <- get_legend(bar_plot)
    bar_plot <- bar_plot + theme(legend.position='none')
    pdf(file=paste0('barplot_site_groups_normalized_',i,'_legend.pdf', sep = ""), width = 22, height = 17)
    print(as_ggplot(bar_plot_legend))
    dev.off()
    pdf(file=paste0('barplot_site_groups_normalized_',i,'.pdf', sep = ""), width = 11, height = 8.5)
    print(bar_plot)
    dev.off()
  }
  
} 

if (replicateFlag == TRUE) {
  setwd(paste0(as.character(args[1]), "/nonnormalized", sep = ""))
  phylo_combined_reps <- merge_samples(phylo_combined, "replicates")
  for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    bar_plot <- plot_bar(phylo_combined_reps, fill = i) +
      geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack")
    bar_plot_legend <- get_legend(bar_plot)
    bar_plot <- bar_plot + theme(legend.position='none')
    pdf(file=paste0('barplot_replicate_groups_nonnormalized_',i,'_legend.pdf', sep = ""), width = 22, height = 17)
    print(as_ggplot(bar_plot_legend))
    dev.off()
    pdf(file=paste0('barplot_replicate_groups_nonnormalized_',i,'.pdf', sep = ""), width = 11, height = 8.5)
    print(bar_plot)
    dev.off()
  }
  setwd(paste0(as.character(args[1]), "/normalized", sep = ""))
  newtotal = median(sample_sums(phylo_combined_reps))
  standf = function(x, t=total) round(t * (x / sum(x)))
  phylo_combined_norm_reps <- transform_sample_counts(phylo_combined_reps, standf)
  for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    bar_plot <- plot_bar(phylo_combined_norm_reps, fill = i) +
      geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack")
    bar_plot_legend <- get_legend(bar_plot)
    bar_plot <- bar_plot + theme(legend.position='none')
    pdf(file=paste0('barplot_replicate_groups_normalized_',i,'_legend.pdf', sep = ""), width = 22, height = 17)
    print(as_ggplot(bar_plot_legend))
    dev.off()
    pdf(file=paste0('barplot_replicate_groups_normalized_',i,'.pdf', sep = ""), width = 11, height = 8.5)
    print(bar_plot)
    dev.off()
  }
  
}


#Heatmaps
setwd(paste0(as.character(args[1]), "/nonnormalized", sep = ""))
heatmap_all <- plot_heatmap(phylo_combined, method = "NMDS", distance = "bray", low = "beige", high = "red", na.value = "beige")
heatmap_combined <- plot_heatmap(phylo_combined_reps, method = "NMDS", distance = "bray", low = "beige", high = "red", na.value = "beige")
pdf(file='heatmap_nonnormalized.pdf', width = 11, height = 8.5)
heatmap_all
dev.off()
pdf(file='heatmap_grouped_nonnormalized.pdf', width = 11, height = 8.5)
heatmap_combined
dev.off()
#GO BACK and put in the if statement section (since it will overwrite if both sites and replicates are provided)

setwd(paste0(as.character(args[1]), "/normalized", sep = ""))
heatmap_all <- plot_heatmap(phylo_combined_norm, method = "NMDS", distance = "bray", low = "beige", high = "red", na.value = "beige")
heatmap_combined <- plot_heatmap(phylo_combined_norm_reps, method = "NMDS", distance = "bray", low = "beige", high = "red", na.value = "beige")
pdf(file='heatmap_normalized.pdf', width = 11, height = 8.5)
heatmap_all
dev.off()
pdf(file='heatmap_grouped_normalized.pdf', width = 11, height = 8.5)
heatmap_combined
dev.off()

#heatmap abundance filtered
# phylo_combined_abund <- filter_taxa(phylo_combined, function(x) sum(x > total*(filter_percent/100)) > 0, TRUE)
# phylo_combined_reps_abund <- filter_taxa(phylo_combined_reps, function(x) sum(x > total*(filter_percent/100)) > 0, TRUE)
# phylo_combined_norm_abund <- blah
# phylo_combined_norm_reps_abund <- blah
#   

#Alpha diversity metrics
#from normalized only
setwd(paste0(as.character(args[1]), "/normalized", sep = ""))
alpha <- plot_richness(phylo_combined_norm, measures = c("Observed", "Shannon"))
pdf(file='alphaDiversity_normalized.pdf', width = 11, height = 8.5)
alpha
dev.off()

alpha_group <- plot_richness(phylo_combined_norm, measures = c("Observed", "Shannon"), x = "replicates", color = "groupB") + geom_point(size=4)
pdf(file='alphaDiversity_grouped_normalized.pdf', width = 11, height = 8.5)
alpha_group
dev.off()


#NMDS ordination
#currently nonnormalized
setwd(paste0(as.character(args[1]), "/nonnormalized", sep = ""))

phylo_combined.ord <- ordinate(phylo_combined, "MDS", "bray")
#phylo_combined_reps.ord <- ordinate(phylo_combined_reps, "MDS", "bray")

ord_ASV <- plot_ordination(phylo_combined, phylo_combined.ord, type="taxa", color="Phylum", 
                title="ASVs")
pdf(file='NMDS_ASVs_nonnormalized.pdf', width = 11, height = 8.5)
ord_ASV
dev.off()

# df <- as.data.frame(phylo_combined.ord$vectors)
# sample_metadata <- tibble::rownames_to_column(sample_metadata)
# df <- tibble::rownames_to_column(df)
# df <- left_join(select(sample_metadata, rowname, replicates), df, by = "rowname" )
# find_hull <- function(df) df[chull(df$Axis.1, df$Axis.2), ]
# hulls <- ddply(df, "replicates", find_hull)

ord_sample <- plot_ordination(phylo_combined, phylo_combined.ord, type="samples", color="replicates", title="Samples") + 
  geom_point(size=3) + 
  geom_encircle(aes(fill = replicates), alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  #geom_polygon(data = hulls, aes(x = Axis.1, y = Axis.2, fill=replicates), alpha = 0.5)
pdf(file='NMDS_sample_nonnormalized.pdf', width = 11, height = 8.5)
ord_sample
dev.off()

# plot_ordination(phylo_combined_reps, phylo_combined_reps.ord, type="taxa", color="Phylum", 
#                 title="ASVs")
# plot_ordination(phylo_combined_reps, phylo_combined_reps.ord, type="samples", 
#                 title="Samples") + geom_point(size=3)






