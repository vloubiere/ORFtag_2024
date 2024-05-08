setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(data.table)

meta <- readRDS("Rdata/processed_metadata_ORFtrap.rds")

stats <- meta[file.exists(bam_stats), {
  # Mapping stats
  .c <- fread(bam_stats, fill= T, sep= "\t")
  mapped <- .c[V1=="reads mapped:"][[2]]
  unmapped <- .c[V1=="reads unmapped:"][[2]]
  #Insertions
  .c <- fread(counts_same_strand)
  insertions <- nrow(.c)
  assigned_insertions <- nrow(.c[dist<2e5])
  # Return
  .(`Aligned reads`= formatC(mapped, big.mark = ",", format = "d"),
    `Aligned %`= round(mapped/(mapped+unmapped)*100, 1),
    `Unique insertions`= formatC(insertions, big.mark = ",", format= "d"),
    `Assigned insertions`= formatC(assigned_insertions, big.mark = ",", format= "d"),
    `Assigned %`= round(assigned_insertions/insertions*100, 1))
}, .(sampleID, bam_stats, counts_same_strand)]
stats <- stats[order(as.numeric(gsub(",", "", `Unique insertions`)))]

pdf("pdf/alignment_statistics.pdf", height = 23, width = 10)
gridExtra::grid.table(stats[, -c(2, 3), with= F])
dev.off()