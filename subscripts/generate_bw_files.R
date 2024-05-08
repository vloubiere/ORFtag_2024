setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(data.table)
require(rtracklayer)
source("git_orftrap_1/functions/ORFTRAP_call_hits.R")

meta <- readRDS("Rdata/selected_screens_processed_metadata_ORFtrap.rds")
meta <- unique(meta[, .(species, screen, condition, bed_file)])

# Generate output file for pos (ps) and neg strands (ns) ----
meta[, bw_file_ps:= ifelse(condition=="input", paste0(species, "_input"), paste0(screen, "_", condition)), species]
meta[, bw_file_ps:= paste0("/groups/stark/vloubiere/projects/ORFTRAP_1/db/bw/", bw_file_ps, "_ps.bw")]
meta[, bw_file_ns:= gsub("_ps.bw$", "_ns.bw", bw_file_ps)]

# Merged bw files selected screens ----
meta[, {
  if(any(!file.exists(c(bw_file_ps, bw_file_ns))))
  {
    .c <- lapply(unique(bed_file), import)
    .c <- do.call("c", .c)
    ps <- .c[strand(.c)=="+"]
    ns <- .c[strand(.c)=="-"]
    cov <- coverage(ps)
    rtracklayer::export.bw(cov, bw_file_ps)
    cov <- coverage(ns)
    export.bw(cov, bw_file_ns)
  }
  .SD
}, .(bw_file_ps, bw_file_ns)]
