#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

############----------------------------------------------------##############
############ Compute overlaps RNA-Seq bam file ##############
############----------------------------------------------------##############

# test if there is at least 2 args: if not, return an error
if (length(args)!=4) {
  stop("Please specify:\n
       [required] 1/ Bed file containing collapsed reads\n
       [required] 2/ gtf containing exons with their frame\n
       [required] 3/ gtf containing transcripts and upstream regions\n
       [required] 4/ output .rds file\n")
}

.libPaths(c("/users/vincent.loubiere/R/x86_64-pc-linux-gnu-library/3.6",  "/software/2020/software/r/3.6.2-foss-2018b/lib64/R/library"))
require(Rsamtools)
require(rtracklayer)
require(GenomicRanges)
require(data.table)

# # For tests
# bed <- "db/bed/Activator_sort_rep1_RNA.bed"
# gtfExons <- "db/gtf/exons_phase_mm10.gtf"
# gtfTranscriptsUpstream <- "db/gtf/transcripts_and_upstream_region_mm10.gtf"
# output <- "db/exon_assignment/PTGR_sort-5_rep1_RNA.rds"

# Variables ----
bed <- args[1]
gtfExons <- args[2]
gtfTranscriptsUpstream <- args[3]
output <- args[4]

# Import data ----
dat <- fread(bed,
             col.names = c("seqnames", "start", "end", "frame", "count", "strand"))
dat[, readID:= .I]
exons <- rtracklayer::import(gtfExons)
exons <- as.data.table(exons)
exons <- exons[, .(seqnames, start, end, gene_id, exon_number, strand, phase, exonPhase)]
exons[, exonPhase:= as.integer(exonPhase)]
trup <- rtracklayer::import(gtfTranscriptsUpstream)
trup <- as.data.table(trup)
trup <- trup[, .(seqnames, start, end, gene_id, type, strand)]
setkeyv(dat, c("seqnames", "strand", "start", "end"))
setkeyv(exons, c("seqnames", "strand", "start", "end"))
setkeyv(trup, c("seqnames", "strand", "start", "end"))

# Overlaps with exons
ov <- foverlaps(exons, dat, nomatch= NULL)
ov[, offset:= ifelse(strand=="+", start-i.start, i.end-end)]
ov[, transcriptPhase:= (exonPhase-offset) %% 3]
# ov[offset>0, transcriptPhase:= (exonPhase-offset) %% 3]
unique(ov[offset==-2, .(exonPhase, transcriptPhase)])
unique(ov[offset==-1, .(exonPhase, transcriptPhase)])
unique(ov[offset==1, .(exonPhase, transcriptPhase)])
unique(ov[offset==2, .(exonPhase, transcriptPhase)])

ov[, expectedPhase:= fcase(frame=="f0", 0L,
                           frame=="f1", 2L,
                           frame=="f2", 1L,
                           default= NA)]
ovExon <- ov[, .(seqnames, start, end, frame, count, strand,
                 readID, exon_number, gene_id,
                 phase, exonPhase, transcriptPhase, expectedPhase,
                 type= "exon")]

# Overlaps with other features
ovFeat <- dat[!(readID %in% ov$readID)]
setorderv(trup, "type")
ovFeat[trup, type:= i.type, on= c("seqnames", "start<=end", "end>=start", "strand"), mult= "first"]
ovFeat[is.na(type), type:= "intergenic"]

# Final object
res <- rbind(ovExon, ovFeat, fill= T)
res <- unique(res)
setorderv(res, "readID")
res[, type:= factor(type, c("exon", "transcript", "upstream", "intergenic"))]

# Overlaps with other features
saveRDS(res,
        output)