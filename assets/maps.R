#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
#args[1]<-"/Users/mcallister/Desktop/test_out/test_figs" #FIGURE OUT directory
#args[2]<-"/Users/mcallister/Desktop/chris_test/CP_all_out/sample_metadata_forR.txt" #sample metadata (only pass if lat and long column headers exist)
#args[3]<-TRUE #whether or not there are replicates (and thus identical coordinates)
#args[4]<-FALSE #whether or not there are site labels in sample metadata
#args[5]<-"/Users/mcallister/Desktop/chris_test/CP_all_out/ASV2Taxonomy/CP_all_out_NO_UNKNOWNS_barchart.txt"
#args[6]<-0.5 #Percent filter setting for pies and barcharts
########################################
library("ggplot2")
library("mapdata")
library("dplyr")
library("mapproj")
library("ggrepel")
library("marmap")
library("scatterpie")
library("janitor")

setwd(as.character(args[1]))
theme_set(theme_bw())

sample_metadata <- read.delim(as.character(args[2]), header=TRUE, stringsAsFactors=FALSE)
replicateFlag <- as.logical(args[3])
sitelabelFlag <- as.logical(args[4])

if (sitelabelFlag == TRUE) {
  sites <- select(sample_metadata, long, lat, sites)
  colnames(sites) <- c("long", "lat", "labels")
} else if (replicateFlag == TRUE) {
  sites <- select(sample_metadata, long, lat, replicates)
  colnames(sites) <- c("long", "lat", "labels")
} else {
  sites <- select(sample_metadata, long, lat, Sample)
  colnames(sites) <- c("long", "lat", "labels")
}

sites <- sites %>% group_by(labels) %>% summarise_all(mean)
sites <- sites %>% filter(!is.na(long))
sites <- sites %>% filter(!is.na(lat))

max_long <- max(sites$long)
min_long <- min(sites$long)
max_lat <- max(sites$lat)
min_lat <- min(sites$lat)

world <- map_data("world2Hires")

if (min_long < 0) {
  minlimlong <- min_long+360-4
} else {
  minlimlong <- min_long-4
}

if (max_long < 0) {
  maxlimlong <- max_long+360+4
} else {
  maxlimlong <- max_long+4
}

minlimlat <- min_lat-4
maxlimlat <- max_lat+4


world_filt <- world %>% filter(long>=minlimlong & long<=maxlimlong)
world_filt <- world_filt  %>% filter(lat>=minlimlat & lat<=maxlimlat)
regions_of_interest <- select(world_filt, region)
collapse <- regions_of_interest %>% group_by(region) %>% summarise_all(funs(toString(unique(.))))
region_string <- sapply(collapse, as.character)
world_of_interest <- map_data("world2Hires", region = c(region_string))

### Simple map with data points
m <- ggplot(data = world_of_interest, aes(x=long, y=lat)) +
  geom_polygon(aes(group = group), 
               fill = "cornsilk",
               color = "black") +
  coord_map(xlim = c(min_long-0.5, max_long+0.5),
            ylim = c(min_lat-0.5, max_lat+0.5), 
            clip = "on") +
  geom_point(data = sites, aes(x = long, y = lat), shape = 21, fill = "red") +
  geom_text_repel(data = sites, aes(x = long, y = lat, label = labels), box.padding = 0.5) +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.1), 
        panel.background = element_rect(fill = "aliceblue"))

pdf(file='mapBasic_datapoints.pdf')
m  
dev.off()

### Bathymetric map with data points

if (min_long-0.5 < 0 & max_long+0.5 > 0) {
  antimeridian <- TRUE
} else if (min_long-0.5 > 0 & max_long+0.5 < 0) {
  antimeridian <- TRUE
} else {
  antimeridian <- FALSE
}

if (max_long - min_long > 10) {
  resolution <- 4
} else {
  resolution <- 1
}

bathymap <- getNOAA.bathy(lon1 = min_long-0.5, lon2 = max_long+0.5,
                          lat1 = min_lat-0.5, lat2 = max_lat+0.5,
                          resolution = resolution, keep = TRUE,
                          antimeridian = antimeridian)

m_bath <- autoplot.bathy(bathymap, geom = c("c", "r"), colour = "white", size = 0.1, coast = TRUE) + 
  scale_fill_etopo() +
  geom_point(data = sites, aes(x = long, y = lat), shape = 21, fill = "red") +
  geom_text_repel(data = sites, aes(x = long, y = lat, label = labels), box.padding = 0.5) +
  xlab("Longitude") + ylab("Latitude") + labs(fill = "Elevation (m)") +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.1))

pdf(file='mapBathy_datapoints.pdf')
m_bath
dev.off()

m_bath_sanslegend <- autoplot.bathy(bathymap, geom = c("c", "r"), colour = "white", size = 0.1, coast = TRUE) + 
  scale_fill_etopo() +
  geom_point(data = sites, aes(x = long, y = lat), shape = 21, fill = "red") +
  geom_text_repel(data = sites, aes(x = long, y = lat, label = labels), box.padding = 0.5) +
  xlab("Longitude") + ylab("Latitude") + labs(fill = "Elevation (m)") +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.1),
        legend.position = "none")

pdf(file='mapBathy_sanslegend_datapoints.pdf')
m_bath_sanslegend
dev.off()

#############Simple map with data points and pie charts
filter_percent <- as.numeric(args[6])

abundance_dat <- read.delim(as.character(args[5]), header=TRUE, stringsAsFactors=FALSE, check.names = FALSE)
abundance_dat <- abundance_dat %>% remove_empty("cols")
abundance_dat$readSum <- rowSums(abundance_dat[, -1])
abundance_dat <- abundance_dat %>% mutate_at(vars(-Sample, -readSum), funs(100*./readSum))
abundance_dat <- abundance_dat %>% select(-readSum)

if (sitelabelFlag == TRUE) {  ##STOP HERE
  abundance_meta_join <- left_join(select(sample_metadata, Sample, sites), abundance_dat, by = "Sample")
  abundance_meta_join <- abundance_meta_join %>% select(-Sample)
  rep_group <- abundance_meta_join %>% group_by(sites) %>% 
    summarise_all(funs(mean)) %>% 
    mutate(zzOther = rowSums(select_if(., ~is.numeric(.) & max(.) < filter_percent))) %>%
    select_if(~is.numeric(.) & max(.) >= filter_percent | is.character(.))
  rep_group <- left_join(select(sample_metadata, sites, lat, long), rep_group, by = "sites")
  rep_group <- rep_group %>% filter(!is.na(long))
  rep_group <- rep_group %>% filter(!is.na(lat))
  rep_group <- rep_group %>% group_by(sites) %>%
    summarise_all(funs(as.numeric(unique(.))))
  number_observations <- ncol(rep_group)
  number_observations <- number_observations - 3
} else if (replicateFlag == TRUE) {
  abundance_meta_join <- left_join(select(sample_metadata, Sample, replicates), abundance_dat, by = "Sample")
  abundance_meta_join <- abundance_meta_join %>% select(-Sample)
  rep_group <- abundance_meta_join %>% group_by(replicates) %>% 
    summarise_all(funs(mean)) %>% 
    mutate(zzOther = rowSums(select_if(., ~is.numeric(.) & max(.) < filter_percent))) %>%
    select_if(~is.numeric(.) & max(.) >= filter_percent | is.character(.))
  rep_group <- left_join(select(sample_metadata, replicates, lat, long), rep_group, by = "replicates")
  rep_group <- rep_group %>% filter(!is.na(long))
  rep_group <- rep_group %>% filter(!is.na(lat))
  rep_group <- rep_group %>% group_by(replicates) %>%
    summarise_all(funs(as.numeric(unique(.))))
  number_observations <- ncol(rep_group)
  number_observations <- number_observations - 3
} else {
  abundance_dat_mod <- abundance_dat %>%
    mutate(zzOther = rowSums(select_if(., ~is.numeric(.) & max(.) < filter_percent))) %>%
    select_if(~is.numeric(.) & max(.) >= filter_percent | is.character(.))
  rep_group <- left_join(select(sample_metadata, Sample, lat, long), abundance_dat_mod, by = "Sample")
  rep_group <- rep_group %>% filter(!is.na(long))
  rep_group <- rep_group %>% filter(!is.na(lat))
  number_observations <- ncol(rep_group)
  number_observations <- number_observations - 3
}

colnames(rep_group)[1] <- "labels"

min_long_dim <- max_long - min_long
min_lat_dim <- max_lat - min_lat
dims <- c(min_long_dim, min_lat_dim)
min_dim <- min(dims)
max_dim <- max(dims)

long_mid <- max_long - (min_long_dim/2)
lat_mid <- max_lat - (min_lat_dim/2)

pie_scale <- (min_dim + 0.5) * 5 

plasma_pal <- c(viridis::plasma(n = number_observations - 1), "lightgrey")

legend_title <- paste("Taxonomic Group (>", filter_percent, "%)", sep = "")

m_pies <- ggplot(data = world_of_interest, aes(x=long, y=lat)) +
  geom_polygon(aes(group = group), 
               fill = "cornsilk",
               color = "black") +
  coord_map(xlim = c(min_long-0.5, max_long+0.5),
            ylim = c(min_lat-0.5, max_lat+0.5), 
            clip = "on", projection = "orthographic") +
  geom_text_repel(data = sites, aes(x = long, y = lat, label = labels), box.padding = pie_scale/5, segment.colour = "grey") +
  geom_scatterpie(data = rep_group, aes(x=long, y=lat, group=labels),
                  cols = colnames(rep_group[,4:ncol(rep_group)]),
                  pie_scale = pie_scale,
                  sorted_by_radius = TRUE, color = NA) +
  scale_fill_manual(values = plasma_pal, name = legend_title) +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.1), 
        panel.background = element_rect(fill = "aliceblue"))

pdf(file='mapPies.pdf')
m_pies
dev.off()
  











