setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")

# All protein coding genes ----
transcripts <- rtracklayer::import("db/gtf/unique_protein_coding_transcripts_mm10.gtf")
transcripts <- as.data.table(transcripts)
genes <- transcripts[, .(exon_number= max(exon_number)), gene_id]

# Import data ----
meta <- readRDS("Rdata/selected_screens_processed_metadata_ORFtrap.rds")[species=="mouse" & condition=="input"]
meta[, screen:= factor(screen, c("Activator", "Repressor", "PTGR"))]
dat <- meta[, fread(counts_same_strand), .(screen, Cc, counts_same_strand)]
dat <- na.omit(dat[dist<=2e05])

# Input pie chart ----
genes[, status:= fcase(gene_id %in% dat$gene_id, "Trapped genes",
                       exon_number==1, "Intronless",
                       default= "No insertion")]
tot <- round(sum(genes$status=="Trapped genes")/nrow(genes)*100, 1)
genes[, status:= paste0(status, "(", formatC(.N, big.mark = ","), ")"), status]

# Saturation curves ----
# randomize
set.seed(1)
sat <- dat[sample(.N)]
allInputs <- sat[, .(Nread= min(.I), screen= "All", Cc= "tomato"), gene_id]
allInputs[, fract:= seq(.N)/nrow(genes)]
sat[, idx:= rowid(screen)]
sat <- sat[, .(Nread= min(idx)), .(screen, Cc, gene_id)]
sat[, fract:= seq(.N)/nrow(genes), screen]
sat <- rbind(allInputs, sat)
sat[, screen:= factor(screen, c("All", "Activator", "Repressor", "PTGR"))]
setorderv(sat, "screen")

# Plot ----
pdf("pdf/saturation_insertions_protein_coding_genes.pdf", 6, 4)
par(tcl= -0.2,
    mgp= c(2, 0.5, 0),
    mai= rep(0.82, 4),
    las= 1)
# Pie chart captured genes
pie(table(genes$status),
    main= paste0("Trapped genes\n(", tot, "% of ", formatC(nrow(genes), big.mark = ","), " genes)"),
    col= c("black", "lightgrey", adjustcolor("#EA6868", 0.6)))
# Satutration curves
par(mai= c(1.02, 0.82+1.25, 0.82, 0.42+1.25))
plot(NA, 
     xlim= c(0, .85), 
     ylim= c(0, 1), 
     xlab= "Number of total insertions (10^6)",
     ylab= "Fraction of protein coding genes",
     frame= F,
     main= "Insertions saturation")
sat[, {
  lines(Nread/1e6,
        fract, 
        col= Cc,
        lwd= 2,
        main= "Insertions saturation")
}, .(screen, Cc)]
# Add line taggable genes
abline(h= sum(genes[, exon_number>1])/nrow(genes), 
       lty= 2)
# Legend
unique(sat[, .(Cc, screen)])[, {
  legend("bottomright",
         legend = screen,
         col= Cc,
         lty= 1,
         lwd= 2,
         bty= "n",
         cex= 0.8)
}]
dev.off()
