setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")

# Mouse ----
gtf <- rtracklayer::import("db/gtf/gencode.vM25.basic.annotation.gtf.gz")
# Unique gene IDs
transcripts <- gtf[gtf$transcript_type=="protein_coding" & !is.na(gtf$transcript_type)]
transcripts <- as.data.table(transcripts)
transcripts[, gene_id:= gsub("[.*]..*", "\\1", gene_id)]
transcripts[, transcript_id:= gsub("[.*]..*", "\\1", transcript_id)]
# Max exon number
transcripts[, exon_number:= max(exon_number, na.rm = T), transcript_id]
# Min CDS length
transcripts <- merge(transcripts,
                     transcripts[type=="CDS", .(CDS_length= sum(end-start+1)), transcript_id], 
                     by= "transcript_id")
# Save as gtf
transcripts <- transcripts[type=="transcript", 
                           .(seqnames, start, end, width, strand, 
                             transcript_id, gene_id, exon_number, CDS_length, gene_name)]
transcripts <- unique(transcripts)
setorderv(transcripts, c("seqnames", "start", "end"))
rtracklayer::export(transcripts, 
                    "db/gtf/unique_protein_coding_transcripts_mm10.gtf")

# Human ----
gtf <- rtracklayer::import("db/gtf/gencode.v43.basic.annotation.gtf.gz")
# Unique gene IDs
transcripts <- gtf[gtf$transcript_type=="protein_coding" & !is.na(gtf$transcript_type)]
transcripts <- as.data.table(transcripts)
transcripts[, gene_id:= gsub("[.*]..*", "\\1", gene_id)]
transcripts[, transcript_id:= gsub("[.*]..*", "\\1", transcript_id)]
# Max exon number
transcripts[, exon_number:= max(exon_number, na.rm = T), transcript_id]
# Min CDS length
transcripts <- merge(transcripts,
                     transcripts[type=="CDS", .(CDS_length= sum(end-start+1)), transcript_id], 
                     by= "transcript_id")
# Save as gtf
transcripts <- transcripts[type=="transcript", 
                           .(seqnames, start, end, width, strand, 
                             transcript_id, gene_id, exon_number, CDS_length, gene_name)]
transcripts <- unique(transcripts)
setorderv(transcripts, c("seqnames", "start", "end"))
rtracklayer::export(transcripts, 
                    "db/gtf/unique_protein_coding_transcripts_hg38.gtf")
