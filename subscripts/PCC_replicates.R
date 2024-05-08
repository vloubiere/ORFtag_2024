setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(Hmisc)
require(data.table)
require(GenomicRanges)
require(vlfunctions)

# Import data ----
meta <- readRDS("Rdata/selected_screens_processed_metadata_ORFtrap.rds")[species=="mouse"]
meta[, screen:= factor(screen, c("Activator", "Repressor", "PTGR"))]
setorderv(meta, c("screen", "condition"))
meta[screen=="PTGR" & replicate==3, replicate:= 2]

# Retrieve assigned insertions per gene ----
dat <- meta[, fread(counts_same_strand)[, .N, gene_id], .(screen, Cc, condition, replicate, counts_same_strand)]
dat <- na.omit(dat)

mat <- dcast(dat, gene_id~screen+condition+replicate, value.var = "N", fun.aggregate = sum)
mat <- as.matrix(mat, 1)
mat <- log2(mat+1)
col <- dat[unlist(tstrsplit(colnames(mat), "_", keep= 1)), Cc, on="screen", mult= "first"]
cor <- varclus(mat, similarity = "pearson", trans = "none")

# Plot scatterplots to chose best screen ----
pdf("pdf/Correlations_select_screens.pdf", 4, 4)
par(mar = c(5,5,2,2),
    tcl= -0.2,
    mgp= c(2,0.5,0),
    las= 1)
plot(cor, hang = -1, names.arg = NULL, xlab = NA, las = 1, ylab = "Pearson correlation coefficient")
points(1:ncol(mat),
       rep(0,ncol(mat)),
       col = col[cor$hclust$order],
       pch = ifelse(grepl("input", colnames(mat)), 17, 19)[cor$hclust$order])
dat[, {
  reps <- CJ(V1= replicate, V2= replicate, unique = T)
  reps <- reps[V2>V1]
  .c <- .SD
  .c <- reps[, merge(.c[replicate==V1, .(gene_id, x= log2(N+1))],
                     .c[replicate==V2, .(gene_id, y= log2(N+1))],
                     by= "gene_id"), .(V1, V2)]
  .c[, PCC:= cor.test(x, y)$estimate, .(V1, V2)]
  .c[, {
    smoothScatter(x, y, 
                  xlab= paste0("Counts rep. ", V1, " (log2)"),
                  ylab= paste0("Counts rep. ", V2, " (log2)"),
                  main= paste(screen, condition))
    legend("topleft",
           legend= paste0("PCC= ", round(PCC, 2)),
           bty= "n")
    .SD
  }, .(V1, V2, PCC)]
  .SD
}, .(screen, condition)]
dev.off()