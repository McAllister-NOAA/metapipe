# MetaPipe

ALL SUPPORT FOR THIS PIPELINE HAS NOW BEEN TRANSFERRED TO REVAMP (URL).

Metabarcoding Pipeline - DEVELOPMENT REPO

Created by Sean McAllister, Chris Paight, Emily Norton, and Matt Galaska

Developed by the Ocean Molecular Ecology Group at <a class="ui-tooltip" title="Pacific Marine Environmental Laboratory"><span style="cursor: help;">PMEL</span></a> (<a class="ui-tooltip" title="National Oceanic and Atmospheric Administration"><span style="cursor: help;">NOAA</span></a>) in cooperation with <a class="ui-tooltip" title="Cooperative Institute for Climate, Ocean, & Ecosystem Studies"><span style="cursor: help;">CICOES</span></a> (<a class="ui-tooltip" title="University of Washington"><span style="cursor: help;">UW</span></a>)

**MetaPipe employs several programs internally, which should also be cited:**

* dada2
* cutadapt
* blastn
* KRONA
* taxonkit
* phyloseq
* vegan

## Installation
### Easy Installation
Docker TBD

### Installation outside of Docker

```
git clone https://github.com/McAllister-NOAA/metapipe.git

```

Add ```metapipe.sh``` to your PATH

Install the necessary dependencies, stand alone packages and in R (below)
##### Dependencies
* R (Rscript)
* dada2 – v.1.14.1
* dbplyr - v.1.4.2
* vegan - v.2.5-6
* mapping packages
* ggplot2 - v.3.3.0
* blastn - v.x
* cutadapt - v.2.8
* taxonkit - v.0.5.0
* krona
* subtree (taxonomy...)
* perl : List::MoreUtils
* And many more...

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
