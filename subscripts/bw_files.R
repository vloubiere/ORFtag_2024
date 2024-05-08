require(rtracklayer)

meta <- readRDS("Rdata/processed_metadata_ORFtrap.rds")
meta <- meta[!is.na(sel)]
meta[, bw_file_ps:= paste0("db/bw/", sel, "_selected_ps.bw")]
meta[, bw_file_ns:= paste0("db/bw/", sel, "_selected_ns.bw")]
meta[, {
  if(any(!file.exists(c(bw_file_ps, bw_file_ns))))
  {
    .c <- lapply(bed_file, import)
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

