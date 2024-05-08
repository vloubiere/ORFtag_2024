#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

############----------------------------------------------------##############
############ Collapse RNA-Seq bam file ##############
############----------------------------------------------------##############

# test if there is at least 2 args: if not, return an error
if (length(args)!=2) {
  stop("Please specify:\n
       [required] 1/ Aligned bam file in which frame has been appended to readID\n
       [required] 2/ Output bed file\n")
}

.libPaths(c("/users/vincent.loubiere/R/x86_64-pc-linux-gnu-library/3.6",  "/software/2020/software/r/3.6.2-foss-2018b/lib64/R/library"))
require(Rsamtools)
require(rtracklayer)
require(GenomicRanges)
require(data.table)

# For tests
# bam <- "db/bam/PTGR_sort-5_rep1_RNA.bam"
# bed <- "db/bed/PTGR_sort-5_rep1_RNA.bed"

# Variables ----
bam <- args[1]
bed <- args[2]

# Import data ----
.c <- Rsamtools::scanBam(bam,
                         param = Rsamtools::ScanBamParam(what= c("qname", "rname", "strand", "pos", "qwidth", "mapq")))
.c <- .c[[1]]
.c <- as.data.table(.c)
setnames(.c,
         c("rname", "pos"),
         c("seqnames", "start"))
# Format
.c[, end:= start+qwidth-1]
.c[, frame:= tstrsplit(qname, "_", keep= 2)]
setorderv(.c,
          c("seqnames", "start", "end", "strand"))

# Alignment statistics ----
stats <- .c[, .(total= .N,
                aligned= sum(mapq>=30, na.rm= T),
                aligned_unique= nrow(unique(.SD[mapq>=30, .(seqnames, start, end, strand)]))), keyby= frame]
fwrite(stats,
       gsub(".bed$", "_stats.txt", bed),
       sep= "\t",
       na = NA)

# Select clean reads ----
coll <- .c[mapq>=30, .(score= .N), .(seqnames, start, end, name= frame, strand)]
setcolorder(coll,
            names(coll)[c(1,2,3,4,6,5)])
fwrite(coll, 
       bed,
       col.names = F,
       row.names = F,
       quote= F,
       sep= "\t",
       na = NA)