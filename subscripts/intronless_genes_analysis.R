setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(vlfunctions)
require(AnnotationDbi)
require(org.Mm.eg.db)

# Import ----
gtf <- rtracklayer::import("db/gtf/gencode.vM25.basic.annotation.gtf.gz")
coding <- gtf[gtf$transcript_type=="protein_coding" & !is.na(gtf$transcript_type)]

# Mark intronic genes as hits ----
genes <- as.data.table(coding)
genes[, hit:= max(exon_number, na.rm= T)==1, gene_id]
genes <- unique(genes[, .(gene_id= gsub("(.*)[.].*", "\\1", gene_id), gene_name, hit)])

# Retrieve protein domains ----
domains <- select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79, 
                  keys = genes$gene_id, 
                  keytype = "GENEID", 
                  columns = "PROTEINDOMAINID")
domains <- unique(na.omit(as.data.table(domains)))
setnames(domains, c("gene_id",  "ID"))

# Import protein domains, restrict to family and merge----
extra <- fread("db/public_db/protein_domains_description_uniqID.txt")
extra <- extra[type=="Family"]
dat <- merge(domains, 
             extra, 
             by.x= "ID", 
             by.y= "ID")
dat <- unique(dat[, .(gene_id, description)])
dat[, hit:= gene_id %in% genes[(hit), gene_id]]
dat[genes, gene_name:= i.gene_name, on= "gene_id"]

# Enrichment ----
dat[, gene_total:= length(unique(genes$gene_id))]
dat[, hit_total:= length(unique(genes[(hit), gene_id]))]
dat <- dat[, .(hit_counts= sum(hit), 
               annot_counts= .N, 
               genes= paste0(unique(gene_name[(hit)]), collapse = ",")), 
           .(description, gene_total, hit_total)]
dat[, c("estimate", "pval"):= {
  .t <- matrix(c(hit_counts,
                 hit_total-hit_counts,
                 annot_counts-hit_counts,
                 gene_total-(hit_total+annot_counts-hit_counts)), 
               ncol= 2)
  fisher.test(.t, alternative = "greater")[c("estimate", "p.value")]
}, .(gene_total, hit_total, hit_counts, annot_counts)]
dat[, padj:= p.adjust(pval, "fdr")]
dat[, log2OR:= log2(estimate)]

# Select top padj for each description and cap log2OR ----
dat <- dat[, .SD[which.min(padj)], description]
dat[, log2OR:= ifelse(log2OR==Inf, max(log2OR[is.finite(log2OR)]), log2OR)]
setorderv(dat, c("log2OR", "padj"), c(-1,1))
dat$pval <- dat$estimate <- NULL
dat <- dat[padj<0.05 & hit_counts>=3 & annot_counts>5]

# Plot enrichments individual screens ----
dat[, padj:= ifelse(padj==0, min(padj[padj>0]), padj)]
setorderv(dat, "log2OR", -1)
dat[, description:= paste0(description, " (", hit_counts, "/", annot_counts, ")")]
breaks <- range(-log10(dat$padj))
if(length(unique(breaks))==1)
  breaks <- breaks+c(-0.5, 0.5)
Cc <- circlize::colorRamp2(breaks, c("blue", "red"))

# Plot ----
pdf("pdf/Intronless_genes_protein_domains_enrichment.pdf", 
    width = 8,
    height =  7)
par(las= 2,
    tcl= -0.2,
    mgp= c(2, 0.25, 0),
    mar= c(25,15,4,0.5))
bar <- barplot(dat$log2OR,
               col= Cc(-log10(dat$padj)),
               border= NA,
               ylab= "Odd ratio (log2)", 
               xaxt= "n")
vl_tilt_xaxis(bar, labels= dat$description)
vl_heatkey(breaks, 
           col= c("blue", "red"),
           left= par("usr")[2]-strwidth("M", cex= 8), 
           top = par("usr")[4]+strwidth("M", cex= 2), 
           main = "padj (-log10)")
dev.off()