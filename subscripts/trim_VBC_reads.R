setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(data.table)
require(parallel)
require(vlfunctions)

# Import metadata ----
meta <- fread("Rdata/RNA_metadata_v2.txt")

# Parse Metadata ----
if(!all(meta[, sampleID==paste0(screen, "_", condition, "_rep", replicate)]))
  stop("sampleID should be the catenation of screen, condtion and replicate, so that re-sequenced samples are correctly collapsed")

# Generate output paths ----
meta[, fq1:= paste0("db/fastq/", gsub(".bam", "", basename(bam_path)), "_", make.unique(sampleID), "_1.fq.gz")]
meta[, fq1_frame:= gsub(".fq.gz$", "_trimmed.fq.gz", fq1)]
meta[, fq2:= as.character(NA)] # data in single end
saveRDS(meta,
        "Rdata/RNA_trimmed_metadata.rds")

# Job parameters and modules ----
cores <- 8
mem <- 32

# Load modules
meta[, load_cmd:= paste(c("cd /groups/stark/vloubiere/projects/ORFTRAP_1/",
                          "module load build-env/2020",
                          "module load trim_galore/0.6.0-foss-2018b-python-2.7.15",
                          "module load samtools/1.9-foss-2018b",
                          "module load bowtie2/2.3.4.2-foss-2018b",
                          "module load r/3.6.2-foss-2018b"), collapse = ";")]

# Demultiplex VBC bam file ----
meta[, demultiplex_cmd:= {
  if(!file.exists(fq1))
  {
    BC <- paste0("'^BC:Z:", barcodes, "'")
    fq_prefix <- normalizePath(gsub("_1.fq.gz$", "", fq1), mustWork = F)
    .f <- "git_orftrap_1/functions/demultiplex_se_14.pl"
    cmd <- paste("samtools view -@", cores-1, bam_path, 
                 "| perl", normalizePath(.f), BC, 
                 fq_prefix, "; gzip", paste0(fq_prefix, "_1.fq"))
    if(!is.na(fq2)) # Not used if single end
      cmd <- paste0(cmd, "; gzip ", fq_prefix, "_2.fq")
    cmd
  }
}, .(bam_path, layout, barcodes, fq1, fq2)]

# Trim reads and append frame to their names ----
meta[, trim_cmd:= {
  if(!file.exists(fq1_frame))
  {
    tmp <- gsub(".fq.gz$", ".fq", fq1_frame)
    paste0("zcat ", fq1,
           " | perl git_orftrap_1/functions/process_ORFtag_RNASeq_reads.pl > ",
           tmp, "; gzip ", tmp)
  }
}, .(fq1, fq1_frame)]

# Submit ----
sub <- meta[, {
  cmd <- unique(na.omit(unlist(.SD)))
  cmd <- paste0(cmd, collapse = "; ")
  if(cmd!=unique(load_cmd))
    cmd
}, sampleID, .SDcols= patterns("_cmd$")]
if(nrow(sub))
{
  sub[, {
    vl_bsub(cmd, 
            cores= 8, 
            m = 20, 
            name = "vlloub", 
            t = '1-00:00:00',
            o= "/groups/stark/vloubiere/projects/ORFTRAP_1/db/logs/",
            e= "/groups/stark/vloubiere/projects/ORFTRAP_1/db/logs/")
  }, .(cmd= V1)]
}