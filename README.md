# metapipe
Metabarcoding Pipeline

Created by Sean McAllister, Matt Galaska, Chris Paight

##Required input files
*Raw reads
*Program settings
primer sequence F
primer sequence R
force merge? T/F [F]
DESeq rarefaction T/F [F]
expected insert size (bp)
sequence read length (bp)
location of sample metadata file
taxa depth (for figures)
control samples (positive/negative)
contaminant taxa list
taxid of interest list

*Sample metadata
sample order
sample groups
sample lat/long
replicate indication
chemistry

##Dependencies
R (Rscript)
dada2 – v.1.14.1
dbplyr - v.1.4.2
vegan - v.2.5-6
mapping packages
ggplot2 - v.3.3.0

blastn - v.x
cutadapt - v.2.8
taxonkit - v.0.5.0
krona
subtree (taxonomy...)

perl : List::MoreUtils
...

acknowledge citations for dada2 (including blog post), cutadapt


##Test dataset

#### Legal Disclaimer
*This repository is a software product and is not official communication
of the National Oceanic and Atmospheric Administration (NOAA), or the
United States Department of Commerce (DOC).  All NOAA GitHub project
code is provided on an 'as is' basis and the user assumes responsibility
for its use.  Any claims against the DOC or DOC bureaus stemming from
the use of this GitHub project will be governed by all applicable Federal
law.  Any reference to specific commercial products, processes, or services
by service mark, trademark, manufacturer, or otherwise, does not constitute
or imply their endorsement, recommendation, or favoring by the DOC.
The DOC seal and logo, or the seal and logo of a DOC bureau, shall not
be used in any manner to imply endorsement of any commercial product
or activity by the DOC or the United States Government.*
