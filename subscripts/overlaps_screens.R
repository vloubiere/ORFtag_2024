setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(data.table)
require(vlfunctions)
require(GenomicRanges)

# Import data ----
dat <- readRDS("Rdata/final_table_mouse.rds")[(hit)]

# Plot ----
pdf("pdf/hits_overlaps.pdf")
plot.new()
vn <- VennDetail::venndetail(split(dat$gene_name, dat$screen))
plot(vn, type= "venn", mycol= unique(dat$Cc))
vl_par(mar= c(20,20,6,7))
vl_upset_plot(split(dat$gene_name, dat$screen), grid.hex = 1.5)$inter
dev.off()