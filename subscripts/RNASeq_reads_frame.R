setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(vlfunctions)

# Import PTGR hits ----
PTGR <- fread("db/FC_tables/selected_screens/PTGR_vs_allMergedInputs.txt")[(hit)]
Activator <- fread("db/FC_tables/selected_screens/Activator_vs_allMergedInputs.txt")[(hit)]

# Import RNA-Seq reads ----
meta <- readRDS("Rdata/RNA_processed_metadata.rds")
dat <- meta[, {
  .c <- readRDS(exon_assignment) 
  .c[, end:= as.integer(end)]
  .c
}, .(user, sampleID, screen)]
dat <- dat[sampleID %in% c("PTGR_input_rep1", "PTGR_sort-5_rep1", "Activator_input_rep1", "Activator_sort_rep1")
           & frame!="ambiguous"
           & count>=10]
dat[type=="transcript", type:= "intron"]
dat[, type:= droplevels(type)]
dat[, cdition:= tstrsplit(sampleID, "_", keep= 2)]
dat[, observedPhase:= exonPhase]
# dat[, observedPhase:= transcriptPhase]
dat[screen=="PTGR", hit:= gene_id %in% PTGR$gene_id]
dat[screen=="Activator", hit:= gene_id %in% Activator$gene_id]

# For each unique RNA-Seq reass, check if one in frame transcript exists ----
# If not, return first line
dat <- dat[order(expectedPhase!=observedPhase)]
dat <- dat[, .SD[1], .(sampleID, readID)]

pdf("pdf/review_RNASeq_frame.pdf", 12, height = 7/3*2)
vl_par(mfrow= c(2,5),
       las= 2,
       mai= rep(.6, 4),
       lwd= 2,
       font.main= 1)
dat[, {
  perc <- .SD[, {
    # Pie chart reads features
    vl_pie(type,
           main= paste0(sampleID, " ov. feat.\n(all integrations uniq)"))
    vl_pie(exon_number,
           main= "Exon number")
    # Number of in frame transcripts
    sum(exonPhase==expectedPhase, na.rm= T)/sum(!is.na(expectedPhase))
  }, sampleID]
  bar <- barplot(perc$V1,
                 ylab= "In-frame transcripts (%)",
                 ylim= c(0,1))
  vl_tilt_xaxis(bar, labels= perc$sampleID)
  
  perc <- .SD[(hit), {
    # Pie chart reads features
    vl_pie(type,
           main= paste0(sampleID, " ov. feat.\n(hit genes uniq)"))
    vl_pie(exon_number,
           main= "Exon number")
    # Number of in frame transcripts
    sum(exonPhase==expectedPhase, na.rm= T)/sum(!is.na(expectedPhase))
  }, sampleID]
  bar <- barplot(perc$V1,
                 ylab= "In-frame transcripts (%)",
                 ylim= c(0,1))
  vl_tilt_xaxis(bar, labels= perc$sampleID)
}, user]
# Same but using sum read counts instead of uniq reads
dat[, {
  perc <- .SD[, {
    # Pie chart reads features
    vl_pie(rep(type, count),
           main= paste0(sampleID, " ov. feat.\n(all integrations total)"))
    vl_pie(exon_number,
           main= "Exon number")
    # Number of in frame transcripts
    sum(rep(exonPhase==expectedPhase, count), na.rm= T)/sum(!is.na(rep(expectedPhase, count)))
  }, sampleID]
  bar <- barplot(perc$V1,
                 ylab= "In-frame transcripts (%)",
                 ylim= c(0,1))
  vl_tilt_xaxis(bar, labels= perc$sampleID)
  
  perc <- .SD[(hit), {
    # Pie chart reads features
    vl_pie(rep(type, count),
           main= paste0(sampleID, " ov. feat.\n(hit genes total)"))
    vl_pie(exon_number,
           main= "Exon number")
    # Number of in frame transcripts
    sum(rep(exonPhase==expectedPhase, count), na.rm= T)/sum(!is.na(rep(expectedPhase, count)))
  }, sampleID]
  bar <- barplot(perc$V1,
                 ylab= "In-frame transcripts (%)",
                 ylim= c(0,1))
  vl_tilt_xaxis(bar, labels= perc$sampleID)
}, user]
dev.off()