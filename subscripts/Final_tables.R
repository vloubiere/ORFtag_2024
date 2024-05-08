setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")

# Import metadata ----
meta <- readRDS("Rdata/selected_screens_processed_metadata_ORFtrap.rds")
meta <- meta[screen!="hRepressor"]
meta <- meta[, .(FC= paste0("db/FC_tables/selected_screens/", paste0(screen, "_vs_allMergedInputs.txt")),
                 rev= paste0("db/FC_tables/selected_screens/", paste0(screen, "_rev_counts_vs_allMergedInputs.txt"))), .(species, screen, Cc)]
dat <- meta[, {
  # Hits
  .c <- fread(FC)
  # Rev enrichment
  rev <- fread(rev)[(hit), gene_id]
  .c[gene_id %in% rev, Comment:= "Enriched for reversed integrations"]
  .c
}, .(species, screen, Cc)]

# Mouse ----
mouse <- dat[species=="mouse", !"species"]
mouse[, screen:= factor(screen, c("Activator", "Repressor", "PTGR"))]
setorderv(mouse, "screen")
# Add validations
sel <- fread("Rdata/validations_selected_hits.txt")
mouse[sel, validation:= i.name, on= c("screen", "gene_name==gene")]
saveRDS(mouse, "Rdata/final_table_mouse.rds")

# Save excel files
xlmm <- copy(mouse)
setnames(xlmm,
         c("padj", "Comment"),
         c("FDR", "comment"))
xlmm$Cc <- xlmm$validation <- NULL
xlmm[, {
  fwrite(.SD,
         paste0("db/excel_tables/", screen, "_screen.txt"),
         sep= "\t",
         quote= F,
         na= "#N/A")
}, screen]

# Human ----
human <- dat[species=="human", !"species"]
human[, screen:= switch(screen, "hActivator"= "Activator", "hRepressor"= "Repressor"), screen]
human[, screen:= factor(screen, c("Activator", "Repressor"))]
human[, screen:= droplevels(screen)]
setorderv(human, "screen")
saveRDS(human, "Rdata/final_table_human.rds")

# Save excel files
xlhs <- copy(human)
setnames(xlhs,
         c("padj", "Comment"),
         c("FDR", "comment"))
xlhs$Cc <- NULL
xlhs[, {
  fwrite(.SD,
         paste0("db/excel_tables/Human_", screen, "_screen.txt"),
         sep= "\t",
         quote= F,
         na= "#N/A")
}, screen]
