setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")

# Import human annotations and all hits ----
dat <- readRDS("Rdata/Annotations_hs_homologs.rds")
mm <- readRDS("Rdata/final_table_mouse.rds")[screen=="Activator"]

# Add classes for overlaps ----
dat[, `p<0.05`:= fcase(mouse_ENSEMBL %in% mm[padj<0.05 & log2OR>1, gene_id] & (Alerasool), "Shared",
                       mouse_ENSEMBL %in% mm[padj<0.05 & log2OR>1, gene_id], "ORFtag specific",
                       (Alerasool), "ORFeome specific",
                       default = "None")]
dat[, `p<0.001`:= fcase(mouse_ENSEMBL %in% mm[padj<0.001 & log2OR>1, gene_id] & (Alerasool), "Shared",
                       mouse_ENSEMBL %in% mm[padj<0.001 & log2OR>1, gene_id], "ORFtag specific",
                       (Alerasool), "ORFeome specific",
                       default = "None")]

# Melt data before plotting ----
.m <- melt(dat,
           c("gene_id", "TF", "AD"),
           c("p<0.05", "p<0.001"),
           value.name = "hit")
.m[, hit:= paste0(hit, " (n= ", formatC(.N, big.mark = ","), ")"), .(variable, hit)]
.m <- melt(.m,
           measure.vars = c("TF", "AD"),
           variable.name = "class")
res <- .m[, {
  .Hit <- hit
  .Class <- value
  .SD[, fisher.test(table(.Hit==hit, .Class))[c("estimate", "p.value")], hit]
}, .(variable, class)]
res <- res[!grepl("None", hit)]
setorderv(res,
          c("variable", "hit"),
          order = c(-1, 1))

# Plots ----
Cc <- c("#776CAE", "#76BD8C", "#CF82B7")
pdf("pdf/review_Alerasool_DelRosso_enrichments.pdf", 5.25, 5)
vl_par(mfrow= c(2,2), lwd= .5)
res[, {
  bar <- barplot(estimate,
                 ylab= "Odd ratio",
                 main= paste0(class, " (", variable, ")"),
                 col= adjustcolor(Cc, .6))
  abline(h= 1,
         lty= "11")
  vl_tilt_xaxis(bar,
                labels = hit)
  text(bar,
       estimate+strheight("M"),
       formatC(p.value, format = "e", digits = 1),
       cex= 6/12,
       offset= 3,
       xpd= T)
}, .(variable, class)]
dev.off()