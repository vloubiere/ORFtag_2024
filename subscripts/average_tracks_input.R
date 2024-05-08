setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")

# All protein coding transcripts ----
transcripts <- rtracklayer::import("db/gtf/unique_protein_coding_transcripts_mm10.gtf")
transcripts <- as.data.table(transcripts)

# Resize to TSSs ----
transcripts[, start:= ifelse(strand=="-", end, start)]
transcripts[, end:= start]

# Compute coverage ----
ps <- vl_bw_average_track(transcripts[strand=="+"], 
                          tracks = "db/bw/mouse_input_ps.bw", 
                          upstream = 10e3,
                          downstream = 10e3,
                          plot = F,
                          col= "#EA6868", 
                          ignore.strand = F)
ns <- vl_bw_average_track(transcripts[strand=="-"], 
                          tracks = "db/bw/mouse_input_ns.bw", 
                          upstream = 10e3,
                          downstream = 10e3,
                          plot = F,
                          col= "#EA6868", 
                          ignore.strand = F)
ns[, score:= -score]
ns[, bin.x:= -bin.x]
cmb <- rbind(ps, ns)
cmb <- cmb[, .(mean= mean(score), se= sd(abs(score))/sqrt(.N)), .(bin.x, name, col)]

# Plot ----
pdf("pdf/average_tracks_input.pdf", 4.5, 3.25)
par(mar= c(5,6,2,2),
    mgp= c(3,0.5,0),
    tcl= -0.2,
    las= 1)
plot(NA,
     xlim= c(-1e4, 1e4),
     ylim= c(-2e-3, 2e-3),
     xlab= "Genomic distance",
     ylab= NA,
     xaxt= "n",
     yaxt= "n",
     frame= F)
title(ylab= "Coverage", line = 3)
axis(1, 
     seq(-10000, 10000, 5000),
     c("-10kb", "-5kb", "TSS", "5kb", "10kb"))
axis(2, 
     seq(-2e-3, 2e-3, 1e-3),
     seq(-2e-3, 2e-3, 1e-3))
cmb[, {
  polygon(c(bin.x, rev(bin.x)), 
          y= c(mean-se, rev(mean+se)) ,
          col= adjustcolor(col, 0.5),
          border= NA)
  lines(bin.x, mean, col= col)
}, .(name, col)]
dev.off()