setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(data.table)

# Import data ----
dat <- readRDS("Rdata/final_table_mouse.rds")
dat[, padj:= ifelse(padj==0, min(padj[padj>0], na.rm= T), padj), screen]
dat[(!hit), Cc:= "lightgrey"]

# Validated candidates ----
dat[!is.na(validation) & validation!=gene_name, Cc:= "lightsalmon"]
dat <- dat[order(screen, Cc!="lightgrey", Cc=="lightsalmon")]
dat[Cc!="lightsalmon", Cc:= adjustcolor(Cc, 0.5)]

# Plot ----
pdf("pdf/review_volcanos_screens_with_numbers_and_revEnrichment.pdf",
    width = 4.2,
    height =  3.75)
par(mar=c (4,4,2,2), 
    las= 1, 
    tcl= -0.2,
    mgp= c(2.25, 0.5, 0))
dat[, {
  plot(log2OR, 
       -log10(padj),
       xlab= "Odd ratio (log2)",
       ylab= "FDR (-log10)",
       col= Cc,
       pch= ifelse(is.na(Comment), 16, 17),
       main= screen,
       xlim= c(-4, switch(as.character(screen), "PTGR"= 6, 8)),
       ylim= c(0, 300),
       frame= F, xpd= T)
  .SD[validation==gene_name, {
    points(log2OR,
           -log10(padj))
  }]
  .SD[!is.na(validation), {
    text(log2OR, 
         -log10(padj), 
         validation, 
         pos= 4, 
         offset= 1,
         xpd= T)
    segments(log2OR+strwidth("M")*0.5, 
             -log10(padj),
             log2OR+strwidth("M")*1.25, 
             -log10(padj))
  }]
  .SD
}, screen]
dev.off()

