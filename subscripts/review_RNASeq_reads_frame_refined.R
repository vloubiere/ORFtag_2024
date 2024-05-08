setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(vlfunctions)

# Import RNA-Seq reads ----
meta <- readRDS("Rdata/RNA_processed_metadata.rds")
meta <- meta[screen=="Activator"]
dat <- meta[, {
  .c <- readRDS(exon_assignment) 
  .c[, end:= as.integer(end)]
  .c
}, .(user, sampleID, screen)]
dat <- dat[frame!="ambiguous" & count>=10]

# Classes ----
dat[type=="transcript", type:= "intron"]
dat[, type:= droplevels(type)]
dat[, cdition:= tstrsplit(sampleID, "_", keep= 2)]
# dat[, observedPhase:= exonPhase]
dat[, observedPhase:= transcriptPhase]

# For each unique RNA-Seq redas, check if one in frame transcript exists ----
# If not, return first line
dat <- dat[order(expectedPhase!=observedPhase)]
dat <- dat[, .SD[1], .(sampleID, readID)]
dat[, exon_number:= as.numeric(exon_number)]
dat[, exon_class:= cut(exon_number, c(-1, 1, 2, 5, Inf), labels= c("1", "2", "3-5", ">5"))]
dat[, exon_class:= factor(exon_class, c("1", "2", "3-5", ">5"))]

# Plot ----
pdf("pdf/review_RNASeq_frame.pdf", 4.8, 7)
vl_par(mfrow= c(3,2),
       las= 2,
       lwd= 0.5,
       font.main= 1)
dat[, {
  par(mai= rep(.8, 4))
  perc <- .SD[, {
    # Pie chart reads features
    vl_pie(rep(type, count),
           main= paste0(sampleID, "\nreads overlap. feat."),
           labels= "p")
    vl_pie(rep(exon_class, count),
           main= paste0(sampleID, " exon number"),
           labels= "p")
    # Number of in frame transcripts
    sum(rep(exonPhase==expectedPhase, count), na.rm= T)/sum(!is.na(rep(expectedPhase, count)))*100
  }, sampleID]
  par(mai= c(.8,1,.8,1))
  bar <- barplot(perc$V1,
                 ylab= "In-frame transcripts (%)",
                 ylim= c(0,100))
  vl_tilt_xaxis(bar,
                labels= perc$sampleID)
  text(bar,
       perc$V1,
       round(perc$V1, 1),
       pos= 3,
       cex= 0.6,
       xpd= T)
}]
dev.off()