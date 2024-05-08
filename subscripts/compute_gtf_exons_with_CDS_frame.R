setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(GenomicRanges)

# Import Mouse exons and CDSs ----
# gtf <- rtracklayer::import("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.basic.annotation.gtf.gz")
gtf <- rtracklayer::import("db/gtf/gencode.vM25.basic.annotation.gtf.gz")
gtf <- as.data.table(gtf)
gtf <- gtf[transcript_type=="protein_coding"]
gtf[, gene_id:= gsub("[.*]..*", "\\1", gene_id)]
gtf[, exon_id:= gsub("[.*]..*", "\\1", exon_id)]
gtf[, transcript_id:= gsub("[.*]..*", "\\1", transcript_id)]
gtf[, exon_number:= as.integer(exon_number)]
gtf[, strand:= as.character(strand)]
gtf[, width:= as.integer(width)]

exons <- gtf[type=="exon"]
exons <- exons[, .(seqnames, start, end, strand, width, exon_number, exon_id, transcript_id, gene_id)]
CDS <- gtf[type=="CDS"]

# Overlap exons and CDSs ----
setorderv(exons, "exon_number")
exons[, overlapsCDS:= exon_id %in% CDS$exon_id]
exons[, startBeforeCDS:= exon_number <= min(exon_number[overlapsCDS]), transcript_id]
exons[CDS, c("offset", "phase", "CDSwidth"):= 
        .(ifelse(strand=="+", i.start-start, end-i.end),
          i.phase,
          i.width), on= "exon_id"]
exons[is.na(CDSwidth), CDSwidth:= 0]
exons[, exonStartMissingAA:= data.table::shift(ceiling(cumsum(CDSwidth)/3), 1, fill = 0),  transcript_id]
exons[, exonStartMissingPerc:= data.table::shift(ceiling(cumsum(CDSwidth)/sum(CDSwidth)*100), 1, fill = 0),  transcript_id]

# Retrieve phase of overlapping CDS and correct for the offset with exon start ----
#         ___  |          ___ |           ___
# phase 0 CAG  | phase 1 CCAG | phase 2 CCCAG
exons[is.na(offset) & (startBeforeCDS), offset:= width]
exons[!is.na(offset), exonPhase:= rev(cumsum(rev(CDSwidth+offset))) %% 3, transcript_id]
exons[transcript_id=="ENSMUST00000000094"]
exons[transcript_id=="ENSMUST00000000028"]
table(exons$phase,
      exons$exonPhase)

# Save ----
exons <- exons[, .(seqnames, start, end, strand,
                   exon_number, exon_id, gene_id,
                   overlapsCDS, startBeforeCDS, exonStartMissingAA, exonStartMissingPerc,
                   exonPhase)]
exons <- unique(exons)
rtracklayer::export(GRanges(exons), 
                    "db/gtf/exons_phase_mm10.gtf")
