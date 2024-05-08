setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")

# Import activator hits ----
act <- fread("db/FC_tables/selected_screens/Activator_vs_allMergedInputs.txt")[(hit)]
actFrame <- list(`Frame+1`= "db/FC_tables/selected_screens/Activator2_frame0_vs_frameSpecificInput.txt",
                 `Frame+2`= "db/FC_tables/selected_screens/Activator2_frame1_vs_frameSpecificInput.txt",
                 `Frame+3`= "db/FC_tables/selected_screens/Activator2_frame2_vs_frameSpecificInput.txt")
actFrame <- lapply(actFrame, fread)
actFrame <- rbindlist(actFrame,
                      idcol = "screen")
actFrame <- actFrame[(hit)]

# Import Frame-specific integrations ----
dat <- readRDS("db/exon_assignment/Activator2_frame_specific_insertions.rds")
dat <- dat[(hit) & dist<2e5]
dat[, screen:= switch(screen,
                      "Frame0"= "Frame+1",
                      "Frame1"= "Frame+2",
                      "Frame2"= "Frame+3",), screen]

# Only keep integrations inside hits and compute stats ----
stats <- dat[, .(correct= sum(compatiblePhase, na.rm = T),
                 total= .N), .(screen, condition)]
stats[, perc:= correct/total*100]
setorderv(stats, "condition")

# Plot ----
pdf("pdf/review_compare_frame_specific_screens_mouse.pdf", 3, 3)
plot.new()
vn <- VennDetail::venndetail(list("3 frames"= act$gene_id,
                                  "Frame specific"= actFrame$gene_id))
plot(vn, type= "venn")
vl_par(mai= c(1.75,1.5,.5, .2))
vl_upset_plot(split(actFrame$gene_id, actFrame$screen),
              cex.grid = 1.2,
              grid.hex = 1.2)
vl_par(mai= c(1,1,.5, .5),
       las= 2)
bar <- barplot(stats$perc,
               col= rev(grey.colors(5)[-c(1, 5)]),
               ylab= "In-frame products (%)")
vl_tilt_xaxis(bar,
              labels = stats[, paste0(screen, " ", condition)],)
dev.off()