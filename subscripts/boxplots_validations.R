setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(vlfunctions)

# Import data ----
dat <- data.table(file= list.files("/groups/stark/nemcko/ORFtrap_paper/db/FCS", ".fcs", full.names = T))
dat[, c("screen", "candidate"):= tstrsplit(basename(file), "_", keep= c(1,2))]
dat[, group:= fcase(candidate %in% c("DsRed", "BFP", "TetR"), "Neg",
                    candidate %in% c("VPR", "Tnrc6bsd", "KRAB"), "Pos",
                    default= "Cand")]
dat[, group:= factor(group, c("Neg", "Pos", "Cand"))]
dat[, Cc:= fcase(group=="Neg", "#706F6F",
                 group=="Pos", "salmon2",
                 screen=="PTGR", "#FDC200",
                 screen=="Activator", "#7EBC87",
                 screen=="Repressor", "#1D71B8")]
dat <- dat[, {
  print(file)
  df <- flowCore::read.FCS(file)
  .c <- as.data.table(df@exprs)
  .c <- .c[`GFP-A`>0, `GFP-A`] 
  .(log2GFP= log2(.c+1))
}, (dat)]

# Sample 25,000 per cdition ----
set.seed(1)
dat <- dat[, {
  if(.N>25000)
    .SD[sample(.N, 25000)] else
      .SD
}, .(screen, candidate, group, Cc)]

# Compute median, pval ----
dat[, med:= median(log2GFP), .(screen, candidate)]
dat <- dat[order(group, med)]
dat[, candidate:= factor(candidate, unique(candidate))]
dat <- dat[, .(log2GFP= .(log2GFP)), .(screen, candidate, group, Cc)]
dat <- merge(dat,
             dat[group=="Neg", .(screen, log2GFP)],
             by= "screen", suffixes= c("", "_neg"))
dat[, pval:= wilcox.test(unlist(log2GFP), unlist(log2GFP_neg), alternative= ifelse(screen=="Activator", "greater", "less"))$p.value, .(screen, candidate)]
dat[, star:= cut(pval, c(-Inf, 1e-3, 1e-2, 5e-1, Inf), c("***", "**", " *", "N.S"))]

# Add n for cditions with <25000 events ----
dat[, N:= lengths(log2GFP)]
dat[N<25000, star:= paste0("n= ", formatC(length(na.omit(unlist(log2GFP))), big.mark = ","), "\n", star), .(candidate, screen)]

# Plots ----
pdf("pdf/validations_fromFCS_new2.pdf", width = 5.3, height =  3.5)
par(tcl= -0.2, 
    mgp= c(2,0.5,0), 
    mar= c(4,8,2,2),
    las= 1,
    lwd= 0.5,
    xpd= T)
dat[, {
  box <- vl_boxplot(log2GFP,
                    horizontal= F,
                    violin= T,
                    whisklty= 1,
                    boxwex = 0.15,
                    viowex = 0.6,
                    tilt.names = T,
                    names= candidate,
                    ylab= "GFP intensity (log2)",
                    viocol= adjustcolor(Cc, 0.5),
                    col= "white",
                    ylim= if(screen=="Activator") NULL else c(7.5, 17))
  text(seq(candidate),
       box$stats[5,],
       star,
       xpd= T,
       pos= 3, 
       cex= ifelse(star=="N.S", 0.5, 0.8))
  text(3:length(candidate),
       par("usr")[3],
       paste0("     ", 1:8),
       pos= 3,
       col= "#E28D8D",
       cex= 0.8)
}, screen]
dev.off()

file.show("pdf/validations_fromFCS_new2.pdf")
