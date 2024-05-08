require(vlfunctions)
require(vlfunctions)

# All mm protein coding genes ----
transcripts <- rtracklayer::import("db/gtf/unique_protein_coding_transcripts_mm10.gtf")
transcripts <- as.data.table(transcripts)
genes <- transcripts[, .(exon_number= max(exon_number), CDS_length= min(as.numeric(CDS_length))), gene_id]
genes <- genes[exon_number>1]

# Retrieve TPMs ----
rnaseq <- read.table(gzfile("/groups/stark/nemcko/work/ORFtrap/general/GSE99971_RNAseq_tpm.txt.gz"), header = T)
rnaseq <- as.data.table(rnaseq)
cols <- c("X18193_mESC_RNAseq_Rep_1.", "X18194_mESC_RNAseq_Rep_2", "X18195_mESC_RNAseq_Rep_3")
rnaseq[, TPM:= rowMeans(.SD), .SDcols= cols]
genes[rnaseq, TPM:= TPM, on= "gene_id==Geneid"]
genes[is.na(TPM), TPM:= 0]
genes[, quant:= cut(TPM, c(-1, 0, quantile(TPM[TPM>0], na.rm = T)[-1]), c("Inactive", paste0("Q", 1:4)))]

# Retrieve trapped genes (at least one count in input) ----
meta <- readRDS("Rdata/selected_screens_processed_metadata_ORFtrap.rds")[species=="mouse" & condition=="input"]
meta[, screen:= factor(screen, c("Activator", "Repressor", "PTGR"))]
dat <- meta[, fread(counts_same_strand), counts_same_strand]
dat <- dat[dist<=2e05 & exon_number<=2] 
genes[, ORFtag:= gene_id %in% dat$gene_id]

# Add hits ----
hits <- readRDS("Rdata/final_table_mouse.rds")[(hit)]
genes[, `ORFtag hits`:= gene_id %in% hits$gene_id]

# Normalized ratios ----
pl <- melt(genes, id.vars = "quant", measure.vars = c("ORFtag", "ORFtag hits"))
pl <- pl[(value), .N, .(variable, quant)]
pl[, ratio:= N/sum(N), variable]
pl[genes[, .(genome= .N/nrow(genes)), quant], norm_ratio:= ratio/i.genome, on= "quant"]
pl[, variable:= factor(variable, c("ORFtag", "ORFtag hits"))]
pl[, quant:= factor(quant, c("Inactive", "Q1", "Q2", "Q3", "Q4"))]
col <- c("#EC6E6C", "#C192C2")

# Plot ----
pdf("pdf/saturation_TPM.pdf",
    width= 4.5,
    height= 4)
par(las= 1,
    tcl= -0.2,
    mgp= c(2,0.5,0))
bar <- barplot(norm_ratio~variable+quant, 
               pl, 
               beside= T,
               border= NA,
               col= adjustcolor(col, 0.6),
               xlab= NA,
               ylab= "Normalized ratio (genome)",
               xaxt= "n")
abline(h= 1, lty= 2)
text(par("usr")[1],
     par("usr")[3]-strheight("M")*2,
     "Genome",
     offset= 0,
     pos= 2, 
     xpd= NA)
rect(bar[1,]-0.5,
     par("usr")[3]-strheight("M")*2.75,
     bar[2,]+0.5,
     par("usr")[3]-strheight("M")*0.75,
     border= NA,
     col= "lightgrey",
     xpd= NA)
text(colMeans(bar),
     par("usr")[3]-strheight("M")*1.75,
     labels = table(genes$quant)[levels(pl$quant)],
     xpd= NA)
axis(1, 
     mean(bar[1:2,1]),
     levels(pl$quant)[1],
     lty= 0,
     line = 1.5)
polygon(c(bar[1,2], bar[2,c(5,5)]),
        par("usr")[3]-strheight("M")*c(4,4,5),
        xpd= T,
        border= NA,
        col= "grey")
text(mean(c(bar[1,2], bar[2,5])),
     par("usr")[3]-strheight("M")*6,
     "Transcriptional level",
     xpd= T)
text(c(bar)-strwidth("M")*0.5,
     par("usr")[3]+strheight("M"),
     pl[, N, keyby= c("quant", "variable")]$N,
     srt= 90,
     pos= 4,
     xpd= NA)
legend("topleft",
       fill= adjustcolor(col, 0.6),
       legend = levels(pl$variable),
       border= NA,
       bty= "n")
dev.off()
