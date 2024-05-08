# ORFtag_2024

This repository contains all custom scripts supporting the conclusions of the study entitled "Proteome-scale tagging and functional screening in mammalian cells by ORFtag".
Authors: Filip Nemčko*, Moritz Himmelsbach*, Vincent Loubiere, Ramesh Yelagandula, Michaela Pagani, Nina Fasching, Julius Brennecke#, Ulrich Elling5#, Alexander Stark#, Stefan L. Ameres#. Nature Methods, 2024.
* These authors contributed equally. # Corresponding authors.

System requirements:
  - Custom scripts generated for this study were written in R (version 4.2.0) using the R studio IDE (https://www.R-project.org/).
  - Custom functions used for ploting and enrichment analysis were wrapped into a R package that was made publicly available at [https://github.com/vloubiere/vlfunction/tree/nature_v2_revised.](https://github.com/vloubiere/vlfunction/tree/ORFtag_2024).
  - No special hardware should be required.

Installation guide:
  - R and RStudio can be downloaded at https://posit.co/download/rstudio-desktop/. Installation time is approximately 20min.
  - The R package containing custome function can be installed using the following command: devtools::install_github("https://github.com/vloubiere/vlfunction/tree/nature_v2_revised"). Installation time is approximately 20min.

Instructions for use:
  - The "main.R" file lists all the scripts that are needed to process the raw sequencing data, analyse the results and plot the panels presented in the manuscript. Script names should be self-explanatory, with additional comments linking them to the corresponding figure panel when relevant.

1/ The "Data processing" section contains the code that was used to process the raw sequencing reads that were deposited on GEO (GSE222193):
  - Align raw RNA-Seq reads, compute gene counts and differential expression analyses.
  - Align CUT&RUN and ChIP-Seq raw sequencing reads, perform peak calling and differential analyses, compute bigwig tracks.
  - Align ATAC-Seq reads, perform peak calling and differential analyses, compute bigwig tracks.

2/ The "Make data tables" section contains the code that was used to generate the tables that were used for all downstream analyses and figures:
  - Tables containing SNV/InDels, SVs and CNVs gDNA alterations. Alignment and calling were performed by Novogene https://www.novogene.com/us-en/ and resultingfiles were made available on GEO.
  - Definition of non-redundant PcG domains in control conditions, using H3K27me3 CUT&RUN.
  - RNA-Seq clustering and resulting table, containing differential expression for all studied conditions. Corresponding data can be found in extended data tables 1 and 2.
  - ATAC-Seq clustering and differential accessibility values for all studied conditions. Corresponding data can be found in extended data table 3.

3/ The "Figures" section contains the code to reproduce the figures from the paper, starting from the data tables generated in the previous section.

4/ The "Supplementary" section contain the code to reproduce the exetended data figures from the paper, starting from the data tables generated in the previous section.

5/ The "Extended data tables" were used to format the extended data tables.

Running the whole set of analyses should be doable within 2 days. Further details regarding the analyses can be found in the Mehods section on the paper.

For any reasonable further request, please contact anne-marie.martinez@igh.cnrs.fr, giacomo.cavalli@igh.cnrs.fr, vincent.loubiere@imp.ac.at
