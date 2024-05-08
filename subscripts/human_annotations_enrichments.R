setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(vlfunctions)

# Import human annotations ----
annot <- readRDS("Rdata/Annotations_hs_homologs.rds")
old <- c("TF", "AD", "Alerasool", "RD", "RBP")
new <- c("TF genes", "Activation domain", "Alerasool et al.", "Repressive domain", "RNA-binding")
setnames(annot, old, new)
annot <- melt(annot,
              id.vars = "mouse_ENSEMBL",
              measure.vars= new,
              variable.name = "class")
annot <- unique(na.omit(annot))
annot[, class:= factor(class, rev(new))]

# Perform enrichment for each screen ----
dat <- readRDS("Rdata/final_table_mouse.rds")
dat[is.na(hit), hit:= F] # We assume that genes with low counts are inactive for the tested function
dat <- dat[, {
  merge(.SD,
        annot,
        by.x= "gene_id",
        by.y= "mouse_ENSEMBL")
}, screen]
res <- dat[, {
  fisher.test(table(value, hit),
              alternative = "greater")[c("estimate", "p.value")]
}, .(class, screen, Cc)]
res[, log2OR:= log2(estimate)]
setorderv(res, "class")
res[log2OR==-Inf, log2OR:= NA]

# Comparative analysis ----
pdf("pdf/mouse_human_features_enrich.pdf",
    width = 2.8,
    height=  3)
par(las= 1,
    mar= c(0.5,14,1,1.5),
    tcl= -0.2,
    mgp= c(2.25, 0.15, 0),
    mfrow= c(3,1),
    oma= c(3,0,0,0))
res[,{
  breaks <- range(-log10(p.value))
  if(length(unique(breaks))==1)
    breaks <- breaks+c(-0.5, 0.5)
  bar <- barplot(-log10(p.value),
                 col= adjustcolor(Cc, 0.6),
                 horiz = T, 
                 names.arg = class,
                 border= NA,
                 xlab= NA,
                 cex.axis = 0.8)
  leg.left <- par("usr")[2]+strwidth("M")*1.5
  segments(-log10(0.01), 
           min(bar)-0.5,
           -log10(0.01),
           max(bar)+0.5,
           lty= 2,
           lwd= 0.5)
  segments(par("usr")[1]-strwidth("M", cex= 11),
           min(bar),
           par("usr")[1]-strwidth("M", cex= 11),
           max(bar),
           xpd= T)
  text(par("usr")[1]-strwidth("M", cex= 12),
       mean(bar),
       screen,
       xpd= T,
       srt= 90)
  print("")
}, .(screen, Cc)]
title(xlab= "Enrichment\np.value (-log10)", xpd= NA)
dev.off()