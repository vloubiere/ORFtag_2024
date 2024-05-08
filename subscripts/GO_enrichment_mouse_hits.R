setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
source("git_orftrap_1/functions/GO_enrichments.R")
require(vlfunctions)
require(AnnotationDbi)
require(rutils)

# Compute GOs enrichment ----
if(!file.exists("db/enrichments/human_annotations_enrichment.txt"))
{
  dat <- readRDS("Rdata/final_table_mouse.rds")[(hit)]
  list_gene_id <- split(dat$gene_id, dat$screen)
  GO_enrichment(list_gene_id = list_gene_id,
                species = "mouse",
                exon_file = "db/gtf/exons_start_mm10.gtf", 
                output_file = "db/enrichments/mouse_annotations_enrichment.txt")
}
dat <- fread("db/enrichments/mouse_annotations_enrichment.txt")
dat[, screen:= factor(screen, c("Activator", "Repressor", "PTGR"))]
setorderv(dat, "screen")


# Plot comparisons between screens ----
dat[, {
  # Plotting parameters
  pdf <- paste0("pdf/", group, "_enrich_compare.pdf")
  pdf(pdf, 
      width = 24.5,
      height= 18)
  par(las= 2,
      mar= c(10,50,4,8),
      mfrow= c(2,2),
      cex= 1)
  .cex <- switch(group, "Localization"= 0.5, "Protein_domains"= 0.3, "GO"= 0.3)
  # Plot for each type
  .SD[, {
    .c <- copy(.SD)
    # Select top 8 enrichments
    .c <- .c[, .SD[description %in% .SD[rowid(screen)<=8, description]]]
    # Order before casting
    .c[, description:= factor(description, levels= unique(description))]
    # dcast
    if(nrow(.c))
    {
      OR <- dcast(.c, description~screen, value.var = "log2OR", drop= T)
      OR <- as.matrix(OR, 1)
      padj <- dcast(.c, description~screen, value.var = "padj", drop= T)
      padj <- as.matrix(padj, 1)
      padj <- -log10(padj)
      vl_balloons_plot(x= OR,
                       color.var= padj,
                       cex.balloons= .cex,
                       x.breaks= seq(0, 8, 2),
                       color.breaks= c(0, 15),
                       col= c("blue", "red"),
                       main= paste0(type, " enrichment"),
                       balloon.size.legend= "OR (log2)",
                       balloon.col.legend= "padj (-log10)", )
    }else
    {
      plot.new()
      title(main= paste0(type))
      text(0.5, 0.5, "No\nEnrichment\nfound")
    }
    
  }, type]
  dev.off()
  pdf
}, group]
