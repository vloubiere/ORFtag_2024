setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(GenomicRanges)

# Mouse ----
# gtf <- rtracklayer::import("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.basic.annotation.gtf.gz")
gtf <- rtracklayer::import("db/gtf/gencode.vM25.basic.annotation.gtf.gz")
gtf <- as.data.table(gtf)
gtf[, gene_id:= gsub("[.*]..*", "\\1", gene_id)]

# Retrieve transcripts ----
transcripts <- gtf[transcript_type=="protein_coding" & type=="transcript"]

# Retrieve upstream regions ----
upstream <- vl_resizeBed(transcripts,
                         "start",
                         upstream = 2e5-1,
                         downstream = 0,
                         genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)

# Make final object
final <- list(transcript= transcripts[, .(seqnames, start, end, strand, gene_id)],
              upstream= upstream[, .(seqnames, start, end, strand, gene_id)])
final <- rbindlist(final,
                   idcol = "type",
                   fill = T)
final[, gene_id:= gsub("[.*]..*", "\\1", gene_id)]
final <- unique(final)
setorderv(final,
          c("seqnames", "start", "end", "strand"))

# Save ----
rtracklayer::export(GRanges(final), 
                    "db/gtf/transcripts_and_upstream_region_mm10.gtf")
