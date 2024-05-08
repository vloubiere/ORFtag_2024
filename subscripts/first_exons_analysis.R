setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")

# Retrieve intronic coding genes ----
if(!exists("gtf"))
  gtf <- rtracklayer::import("db/gtf/gencode.vM25.basic.annotation.gtf.gz")
coding <- gtf[gtf$transcript_type=="protein_coding" & !is.na(gtf$transcript_type)]
coding_genes <- as.data.table(coding)
coding_genes[, intron:= max(exon_number, na.rm= T)>1, transcript_id]
coding_genes <- coding_genes[(intron)]
coding_genes[, exon_id:= gsub("(.*)[.].*", "\\1", exon_id)]

# For each exon, check whether contains CDS ----
exon <- unique(coding_genes[type=="exon", .(exon_id, exon_number, seqnames, start, end, strand)])
# Contains CDS
CDS <- unique(coding_genes[type=="CDS", .(exon_id, cds_start= start, cds_end= end)])
exon <- merge(exon, CDS, by= "exon_id", all.x= T)
exon[, contains_CDS:= !is.na(cds_start)]
# Size of the CDS
exon[, CDS_small:= cds_end-cds_start+1<=(20*3)]

# For each exon, check whether contains PFAM domain ----
# Contain PFAM
pfam <- fread("db/gtf/pfam_ucsc_mm10_2019-09-20.tsv")
# Get the coding coordinates only
pfam <- pfam[, {
  start <- as.numeric(strsplit(chromStarts, ",")[[1]])
  width <- as.numeric(strsplit(blockSizes, ",")[[1]])
  .(dom_width= sum(width),
    start= chromStart+start,
    end= chromStart+start+width)
}, .(seqnames= chrom, chromStart, chromEnd, strand, name)]
# Overlap with exons and compute overlapping width
setkeyv(exon, c("seqnames", "start", "end"))
setkeyv(pfam, c("seqnames", "start", "end"))
ov <- foverlaps(exon, pfam)
ov <- ov[strand==i.strand]
# Compute Overlaps
ov[, ov_start:= apply(.SD, 1, max), .SDcols= c("start", "i.start")]
ov[, ov_end:= apply(.SD, 1, min), .SDcols= c("end", "i.end")]
# For each exon and each domain, compute the fraction of overlap
ov <- ov[, .(dom_width, ov_width= sum(ov_end-ov_start+1)), .(seqnames, chromStart, chromEnd, name, exon_id)]
ov <- ov[ov_width>0.1*dom_width] #Only consider overlapping if >0.25% domain within exon
exon[, contains_domain:= exon_id %in% ov$exon_id]

# Exon status ----
exon[, status:= fcase((!contains_CDS), "No coding sequence",
                      (contains_domain), "Pfam domain-containing",
                      (contains_CDS) & (CDS_small), "Short coding sequence (<=20aa)",
                      (contains_CDS) & !(CDS_small), "Long coding sequence (>20aa)")]
exon[, status:= factor(status, c("No coding sequence", 
                                 "Pfam domain-containing",
                                 "Long coding sequence (>20aa)",
                                 "Short coding sequence (<=20aa)"))]

# Plot exon fraction ----
pl <- table(exon[exon_number==1, status])

# fraction <- pfam[, .(total= .N, exon_1= sum(exon_number==1)), domain]
# fraction[, perc:= exon_1/total*100]
# setorderv(fraction, c("perc", "exon_1"), -1)

pdf("pdf/first_exon_biased_pfam_domains.pdf", 
    height= 5.5)
par(las= 1,
    tcl= -0.2,
    mar= c(5,12,5,12),
    mgp= c(3,0.5,0))
pie(pl, 
    labels = paste0(names(pl), "\n", formatC(pl, big.mark = ",")))
par(mar= c(10,4,10,1),
    mgp= c(2,0.5, 0))
# fraction[perc>50 & total>=10, {
#   bar <- barplot(perc,
#                  ylab= "% of domains within exon 1",
#                  border= NA,
#                  xaxt= "n",
#                  yaxt= "n")
#   vl_tilt_xaxis(bar, labels= domain)
#   axis(2, 
#        seq(0, 100, length.out= 5),
#        seq(0, 100, length.out= 5))
#   text(bar,
#        perc,
#        exon_1,
#        xpd= T,
#        pos= 3,
#        offset= 0.25,
#        cex= 0.8)
# }]
dev.off()
