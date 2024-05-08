setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(data.table)
require(parallel)

CR <- fread("/groups/stark/nemcko/work/cutandrun/mm10/20230621_Gabpa_Kansl3_Pprc1_Zfp574/expFile.txt")
CR <- CR[grepl("parental_aV5.*0IAA|Zfp574.*0IAA", outfile_prefix)]
CR[, bam_path:= "/groups/stark/nemcko/work/cutandrun/mm10/20230621_Gabpa_Kansl3_Pprc1_Zfp574/HNF5HBGXT_1_20230620B_20230621.bam"]
PRO <- fread("/groups/stark/nemcko/work/proseq/mm10/20211114_Hcfc1_Zfp574_3,6,12,24h/experimentFile_with_eBC.txt")
PRO <- PRO[grepl("Zfp574.*0hr|Zfp574.*6hr", outfile_prefix)]
PRO[, bam_path:= "/groups/stark/nemcko/work/proseq/mm10/20211114_Hcfc1_Zfp574_3,6,12,24h/20211114_PROseq_merged.bam"]
meta <- rbind(CR, PRO, fill= T)
meta[, fq1:= paste0("db/geo_submission/review_1/", gsub(".bam$", "", basename(bam_path)), outfile_prefix, "_1.fq.gz")]
meta[, fq2:= paste0("db/geo_submission/review_1/", gsub(".bam$", "", basename(bam_path)), outfile_prefix, "_2.fq.gz")]

# Demultiplex CUT&RUN data
meta[is.na(eBC), demultiplex_cmd:= {
  if(.N!=1)
    stop("Unique bam file should be provided!")
  if(!file.exists(fq1))
  {
    fq_prefix <- normalizePath(gsub("_1.fq.gz$", "", fq1), mustWork = F)
    cmd <- paste("samtools view -@ 5 ", bam_path, 
                 "| perl", 
                 normalizePath("git_orftrap_1/functions/demultiplex_pe_14.pl"),
                 paste0("'^BC:Z:", `#barcode`, "'"), 
                 fq_prefix, 
                 "; gzip", 
                 paste0(fq_prefix, "_1.fq"))
    if(!is.na(fq2))
      cmd <- paste0(cmd, "; gzip ", fq_prefix, "_2.fq")
    cmd
  }
}, .(`#barcode`, fq1, fq2)]

# Demultiplex PROSeq data
meta[!is.na(eBC), demultiplex_cmd:= {
  if(.N!=1)
    stop("Unique bam file should be provided!")
  if(!file.exists(fq1)) # fq2 is not used anyway
  {
    fq_prefix <- normalizePath(gsub("_1.fq.gz$", "", fq1), mustWork = F)
    cmd <- paste("samtools view -@ 5 ", bam_path, 
                 "| perl", 
                 normalizePath("git_orftrap_1/functions/demultiplex_pe_14_PROSeq.pl"),
                 paste0("'^BC:Z:", `#barcode`, "'"),
                 paste0("'^", eBC, "'"),
                 fq_prefix, 
                 "; gzip", 
                 paste0(fq_prefix, "_1.fq"))
    if(!is.na(fq2))
      cmd <- paste0(cmd, "; gzip ", fq_prefix, "_2.fq")
    cmd
  }
}, .(`#barcode`, eBC, fq1, fq2)]

meta[, demultiplex_cmd:= paste0("module load build-env/2020; module load samtools/1.9-foss-2018b; ", demultiplex_cmd)]
meta[, {
  vl_bsub(cmd, 
          cores= 6, 
          m = 32, 
          name = "vlloub", 
          t = '1-00:00:00',
          o= "/groups/stark/vloubiere/projects/ORFTRAP_1/logs/",
          e= "/groups/stark/vloubiere/projects/ORFTRAP_1/logs/")
}, .(cmd= demultiplex_cmd)]