setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(vlfunctions)

# Import insertions data ----
meta <- readRDS("Rdata/selected_screens_processed_metadata_ORFtrap.rds")
dat <- meta[species=="mouse" & condition=="input", {
 .c <- lapply(counts_same_strand, fread)
 .c <- rbindlist(.c)
}, .(screen, condition, replicate)]
dat <- dat[dist<2e5]

# # Add hits column ----
# hits <- readRDS("Rdata/final_table_mouse.rds")
# dat[hits, hit:= i.hit, on= c("gene_id", "screen")]
# dat[hits, gene_name:= i.gene_name, on= c("gene_id", "screen")]
# dat[condition=="input", screen:="input"]

# Import exon IDs and CDS start ----
exons <- rtracklayer::import("db/gtf/exons_phase_mm10.gtf")
exons <- as.data.table(exons)
exons[, startBeforeCDS:= as.logical(startBeforeCDS)]
exons[, exonStartMissingAA:= as.numeric(exonStartMissingAA)]
exons[, exonStartMissingPerc:= as.numeric(exonStartMissingPerc)]

# Match ----
dat[exons, exonStartUpstreamCDS:= ifelse(startBeforeCDS, "Before", "After"), on= "exon_id"]
dat$missingPerc <- exons[dat, min(exonStartMissingPerc), .EACHI, on= "exon_id"]$V1
dat[, missingPerc:= cut(missingPerc, c(-Inf, 0, 10, Inf), c("Full-length protein", ">90% of protein", "<90% of protein"), include.lowest= T)]
setorderv(dat, "missingPerc")
geneCollPerc <- dat[, .SD[1], gene_id]

# Plot ----
pdf("pdf/review_integration_after_initial_codon.pdf", 14, 3.5)
vl_par(mfrow=c(1,2),
       mai= c(.9,2.9,.9,2.9))
dat[, {
  vl_pie(missingPerc,
         main= "Number of integrations")
}]
geneCollPerc[, {
  vl_pie(missingPerc,
         main= "Number of genes")
}]
dev.off()