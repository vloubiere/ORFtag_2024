setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(vlfunctions)

# All mm protein coding genes ----
transcripts <- rtracklayer::import("db/gtf/unique_protein_coding_transcripts_mm10.gtf")
transcripts <- as.data.table(transcripts)
genes <- transcripts[, .(exon_number= max(exon_number), CDS_length= min(as.numeric(CDS_length))), gene_id]
genes <- genes[exon_number>1]
genes[, CDS_length:= fcase(CDS_length>5000, ">5kb",
                           CDS_length>2500, "2.5-5kb",
                           default = "<2.5kb")]

# Retrieve trapped genes (at least one count in input) ----
meta <- readRDS("Rdata/selected_screens_processed_metadata_ORFtrap.rds")[species=="mouse" & condition=="input"]
meta[, screen:= factor(screen, c("Activator", "Repressor", "PTGR"))]
dat <- meta[, fread(counts_same_strand), counts_same_strand]
dat <- dat[dist<=2e05 & exon_number<=2] 
genes[, ORFtag:= gene_id %in% dat$gene_id]

# Add hits ----
hits <- readRDS("Rdata/final_table_mouse.rds")[(hit)]
genes[, `ORFtag hits`:= gene_id %in% hits$gene_id]

# Retrieve human ORFeome ----
hORF <- fread("db/public_db/hORFeome9.1_annotation.tsv.gz")
hORF <- hORF[orf_class=="pcORF", 
             .(gene_id= gsub("(^.*)[.].*$", "\\1", ensembl_gene_id), 
               CDS_size= nchar(cds))] 
hORF[, CDS_length:= fcase(CDS_size>5000, ">5kb",
                          CDS_size>2500, "2.5-5kb",
                          default = "<2.5kb")]
screen <- readxl::read_xlsx("db/public_db/alerasool_2022.xlsx", sheet = 2)
screen <- as.data.table(screen)
ensembl <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                 keys = screen$`Gene ID`,
                                 columns = "ENSEMBL",
                                 keytype = "ENTREZID")
hORF[, trapped:= gene_id %in% ensembl$ENSEMBL]

# Normalized ratios ----
mm <- melt(genes, id.vars = "CDS_length", measure.vars = c("ORFtag", "ORFtag hits"))
mm <- mm[(value), .N, .(variable, CDS_length)]
mm[, ratio:= N/sum(N), variable]
mm[genes[, .(genome= .N/nrow(genes)), CDS_length], norm_ratio:= ratio/i.genome, on= "CDS_length"]

hs <- hORF[(trapped), .N, CDS_length]
hs[, ratio:= N/sum(N)]
hs[hORF[, .(genome= .N/nrow(hORF)), CDS_length], norm_ratio:= ratio/i.genome, on= "CDS_length"]
hs[, variable:= "ORFeome"]

pl <- rbind(mm, hs)
pl[, variable:= factor(variable, c("ORFeome", "ORFtag", "ORFtag hits"))]
pl[, CDS_length:= factor(CDS_length, c("<2.5kb", "2.5-5kb", ">5kb"))]

col <- c("#9B8579", "#EC6E6C", "#C192C2")

# Plot ----
pdf("pdf/CDS_size_genome_ORFeome_ORFtrap_hits_norm_ratio.pdf",
    width= 4.5,
    height= 4)
par(las= 1,
    tcl= -0.2,
    mgp= c(2,0.5,0))
bar <- barplot(norm_ratio~variable+CDS_length, 
               pl, 
               beside= T,
               border= NA,
               col= adjustcolor(col, 0.6),
               xlab= NA,
               ylab= "Normalized ratio (genome)",
               xaxt= "n")
abline(h= 1, lty= 2)
axis(1, bar[2,], levels(pl$CDS_length), lty= 0, line = 1.25)
title(xlab= "ORF length", line = 3)
text(par("usr")[1],
     par("usr")[3]-strheight("M")*2,
     "Genome",
     offset= 0,
     pos= 2, 
     xpd= NA)
text(c(bar)-strwidth("M")*0.5,
     par("usr")[3]+strheight("M"),
     pl[, N, keyby= c("CDS_length", "variable")]$N,
     srt= 90,
     pos= 4,
     xpd= NA, )
rect(bar[1,]-0.5,
     par("usr")[3]-strheight("M")*2.75,
     bar[3,]+0.5,
     par("usr")[3]-strheight("M")*0.75,
     border= NA,
     col= "lightgrey",
     xpd= NA)
text(bar[2,],
     par("usr")[3]-strheight("M")*1.75,
     labels = table(genes$CDS_length)[levels(pl$CDS_length)],
     xpd= NA)
legend("topleft",
       fill= adjustcolor(col, 0.6),
       legend = levels(pl$variable),
       border= NA,
       bty= "n")
dev.off()