# ORFtag_2024

This repository contains all custom scripts supporting the conclusions of the study entitled "Proteome-scale tagging and functional screening in mammalian cells by ORFtag". Nature Methods, 2024.

Authors: Filip Nemčko†, Moritz Himmelsbach†, Vincent Loubiere, Ramesh Yelagandula, Michaela Pagani, Nina Fasching, Julius Brennecke#, Ulrich Elling5#, Alexander Stark#, Stefan L. Ameres#. 
† These authors contributed equally. # Corresponding authors.

System requirements:
  - Custom scripts generated for this study were written in R (version 4.2.0) using the R studio IDE (https://www.R-project.org/).
  - Custom functions used for ploting and enrichment analysis were wrapped into a R package that was made publicly available at [https://github.com/vloubiere/vlfunction/tree/ORFtag_2024.](https://github.com/vloubiere/vlfunction/tree/ORFtag_2024).
  - No special hardware should be required.

Installation guide:
  - R and RStudio can be downloaded at https://posit.co/download/rstudio-desktop/. Installation time is approximately 20min.
  - The R package containing custom functions can be installed using the following command: devtools::install_github("https://github.com/vloubiere/vlfunction/tree/ORFtag_2024"). Installation time is approximately 5min.

Instructions for use:
  - The "main.R" file lists all the scripts that are needed to process the raw sequencing data, analyse the results and plot the panels presented in the manuscript. Script names should be self-explanatory, with additional comments linking them to the corresponding figure panel when relevant.

1/ The "Data processing" section contains the scripts that were used to process the raw NGS sequencing data that was deposited on GEO (GSE225972), corresponding to:
  - Three different functional ORFtag screens in mouse ESCs (activator, repressor, PTGR)
  - Frame-specific activator ORFtag screen
  - ORFtag RNA-Seq reads used to assess the frame of transcripts (activator screen)
  - CUT&RUN data for V5-tagged zfp574
  - PRO-Seq data after rapid depletion of zfp574 (AID cell-line)

2/ The "Figures" section contains the scripts to reproduce the figures from the article, starting from the different data tables generated in the previous section and are also available on GEO (GSE225972) or in the supplementary tables of the article.

3/ The "Supplementary Figures" section contain the scripts to reproduce the supplementary figures, starting from the data tables generated in section 1.

Running the whole set of analyses should be doable within 1 day. Further details regarding the analyses can be found in the Mehods section on the paper.

For any reasonable further request, please contact alexander.stark@imp.ac.at, filip.nemcko@imp.ac.at or vincent.loubiere@imp.ac.at.
