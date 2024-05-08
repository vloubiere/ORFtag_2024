setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
source("git_orftrap_1/functions/GO_enrichments.R")

# Merge all screens
enr <- rbindlist(list(mouse= fread("db/enrichments/mouse_annotations_enrichment.txt")[screen %in% c("Activator", "Repressor")],
                      human= fread("db/enrichments/human_annotations_enrichment.txt")),
                 idcol= "species")
enr[, screen:= paste(species, screen)]
enr[, screen:= factor(screen,
                      c("mouse Activator",
                        "human Activator",
                        "mouse Repressor",
                        "human Repressor"))]

# Plot comparisons between screens ----
files <- enr[, {
  # Plotting parameters
  pdf <- paste0("pdf/human_vs_mouse_hits_", group, "_compare.pdf")
  pdf(pdf, 
      width = 25.5,
      height= 18)
  par(las= 2,
      mar= c(10,50,4,8),
      mfrow= c(2,2),
      cex= 1)
  .cex <- switch(group, "Localization"= 0.5, "Protein_domains"= 0.3, "GO"= 0.3)
  # Plot for each type
  .SD[, {
    .c <- copy(.SD)
    
    # Select TERMs
    .c <- .c[, .SD[description %in% .SD[rowid(screen)<=10, description]]]# Select top 8 enrichments
    
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
                       color_var= padj,
                       cex.balloons= .cex,
                       x_breaks= seq(0, 8, 2),
                       color_breaks= c(0, 15),
                       col= c("blue", "red"),
                       main= paste0(type, " enrichment "),
                       balloon_size_legend= "OR (log2)",
                       balloon_col_legend= "padj (-log10)", )
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
files
