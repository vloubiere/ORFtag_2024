setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(GenomicRanges)

# Mouse ----
# gtf <- rtracklayer::import("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.basic.annotation.gtf.gz")
gtf <- rtracklayer::import("db/gtf/gencode.vM25.basic.annotation.gtf.gz")
exons <- gtf[gtf$transcript_type=="protein_coding" & gtf$type=="exon" & gtf$exon_number>1]
exons <- GenomicRanges::resize(exons, 1, "start")
mcols(exons) <- mcols(exons[, c("gene_id", "gene_name", "mgi_id", "exon_number", "exon_id")])
exons$gene_id <- gsub("[.*]..*", "\\1", exons$gene_id)
exons$exon_id <- gsub("[.*]..*", "\\1", exons$exon_id)
exons <- unique(exons)
rtracklayer::export(exons,
                    "db/gtf/exons_start_mm10.gtf")

# Human ----
# gtf <- rtracklayer::import("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz")
gtf <- rtracklayer::import("db/gtf/gencode.v43.annotation.gtf.gz")
exons <- gtf[gtf$transcript_type=="protein_coding" & gtf$type=="exon" & gtf$exon_number>1]
exons <- GenomicRanges::resize(exons, 1, "start")
mcols(exons) <- mcols(exons[, c("gene_id", "gene_name", "exon_number", "exon_id")])
exons$gene_id <- gsub("[.*]..*", "\\1", exons$gene_id)
exons$exon_id <- gsub("[.*]..*", "\\1", exons$exon_id)
exons <- unique(exons)
rtracklayer::export(exons, 
                    "db/gtf/exons_start_hg38.gtf")
