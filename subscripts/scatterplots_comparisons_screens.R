setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(vlfunctions)

# Import data ----
dat <- readRDS("Rdata/final_table_mouse.rds")
dat[, screen:= factor(screen)]

# Combinations ----
cmb <- CJ(V1= dat$screen, 
          V2= dat$screen,
          unique = T)
cmb <- cmb[as.numeric(V1)<as.numeric(V2)]
cmb <- cmb[,{
  merge(dat[.BY, .(gene_id, padj= -log10(padj), hit, Cc), on= "screen==V1"],
        dat[.BY, .(gene_id, padj= -log10(padj), hit, Cc), on= "screen==V2"],
        by= "gene_id")
}, .(V1, V2)]
cmb <- cmb[(hit.x) | (hit.y)]
cmb[is.na(padj.x), padj.x:= 0]
cmb[is.na(padj.y), padj.y:= 0]
cmb[, padj.x:= ifelse(is.infinite(padj.x),
                      max(padj.x[!is.infinite(padj.x)]),
                      padj.x), .(V1, V2)]
cmb[, padj.y:= ifelse(is.infinite(padj.y),
                      max(padj.y[!is.infinite(padj.y)]),
                      padj.y), .(V1, V2)]
cmb[, Cc:= fcase(hit.x & hit.y, "tomato",
                 hit.x, Cc.x,
                 hit.y, Cc.y)]

# Plot ----
pdf("pdf/scatterplots_comparison_screens_FDR.pdf", 6, 6.5)
layout(matrix(c(1,2,4,3), byrow = T, nrow = 2))
par(mgp= c(2.25,0.5,0),
    las= 1,
    tcl= -0.2)
cmb[, {
  .SD[, {
    plot(padj.x, 
         padj.y,
         col= adjustcolor(Cc, 0.4),
         pch= 16,
         frame= F,
         xlab= paste0(V1, " (-log10FDR)"),
         ylab= paste0(V2, " (-log10FDR)"))
    text(grconvertX(0.6, "npc", "user"),
         grconvertY(0.15, "npc", "user"), 
         sum(Cc==Cc.x),
         col= Cc.x[1])
    text(grconvertX(0.15, "npc", "user"),
         grconvertY(0.15, "npc", "user"), 
         sum(Cc=="tomato"),
         col= "tomato")
    text(grconvertX(0.15, "npc", "user"),
         grconvertY(0.6, "npc", "user"), 
         sum(Cc==Cc.y),
         col= Cc.y[1])
  }]
}, .(V1, V2)]
dev.off()