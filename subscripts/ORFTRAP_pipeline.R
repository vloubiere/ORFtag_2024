setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(data.table)
require(parallel)
require(vlfunctions)

# ORF-Tag pipeline ----
# Takes as input a metadata table consisting of
#   1- Path to VBC bam file
#   2- barcode pairs, separated by "|"
#   3- A sampleID column, cooresponding to a unique combination of screen, condition, replicate (1 combination/sample)
# 1/ Demultiplexes VBC bam according to barcode pairs, and output fastq files 
# 2/ Aligns first reads to the genome of reference 
#    ! -> if multiple bam files have similar sampleID, they are merged at this step!
# 3/ For each bam file, alignment statistics are computed
# 4/ Unique reads with same start and end and mapq>=30 are collapsed using samtools, and saved using bam format
# 5/ For each bam_unique, computes the exact position of unique insertions and save them in a bed file, before computing the closest downstream non-first exon junction (same and reverse strands)
#    ! -> See the bamToBed_and_assign_insertions.R function for further details
#    -> saved as: assignment_table"_same_strand.txt" & assignment_table"_rev_strand.txt"

## Parse Metadata ----
meta <- fread("Rdata/metadata_v5.txt", header = T)
if(!all(meta[, sampleID==paste0(screen, "_", condition, "_rep", replicate)]))
  stop("sampleID should be the catenation of screen, condtion and replicate, so that re-sequenced samples are correctly collapsed")

## Create output folders ----
dir.create("db/", showWarnings = F)
dir.create("db/fastq", showWarnings = F)
dir.create("db/bam", showWarnings = F)
dir.create("db/bam_unique", showWarnings = F)
dir.create("db/bed", showWarnings = F)
dir.create("db/gene_assignment", showWarnings = F)

## Generate output paths ----
meta[, fq1:= paste0("db/fastq/", gsub(".bam", "", basename(bam_path)), "_", make.unique(sampleID), "_1.fq.gz")]
meta[(layout=="PAIRED"), fq2:= gsub("_1.fq.gz$", "_2.fq.gz", fq1)]
meta[, fq1_trim:= gsub(".fq.gz$", "_trimmed.fq.gz", fq1)]
meta[, bam:= paste0("db/bam/", sampleID, ".bam")] # re-sequencing are merged from this step on!
meta[, bam_stats:= gsub(".bam$", "_stats.txt", bam)]
meta[, bam_unique:= paste0("db/bam_unique/", sampleID, "_q30_unique.bam")]
meta[, bed_file:= paste0("db/bed/", sampleID, ".bed")]
meta[, counts_same_strand:= paste0("db/gene_assignment/", sampleID, "_same_strand.txt")]
meta[, counts_rev_strand:= paste0("db/gene_assignment/", sampleID, "_rev_strand.txt")]

## Job parameters and modules ----
cores <- 8
mem <- 32
# Load modules
meta[, load_cmd:= paste(c("module load build-env/2020",
                          "module load trim_galore/0.6.0-foss-2018b-python-2.7.15",
                          "module load samtools/1.9-foss-2018b",
                          "module load bowtie2/2.3.4.2-foss-2018b",
                          "module load r/3.6.2-foss-2018b",
                          paste0("cd ", getwd())), collapse = "; ")]

## Demultiplex VBC bam file ----
meta[, demultiplex_cmd:= {
  if(.N!=1)
    stop("Unique bam file should be provided!")
  if(!file.exists(fq1)) # fq2 is not used anyway
  {
    BC <- paste0("'^BC:Z:", paste0(unlist(tstrsplit(barcodes, "\\|")), collapse = "|^BC:Z:"), "'")
    fq_prefix <- gsub("_1.fq.gz$", "", fq1)
    .f <- fcase(layout=="SINGLE" & sequencer=="NextSeq",
                "/groups/stark/vloubiere/projects/ORFTRAP_1/git_orftrap_1/functions/demultiplex_se_12.pl", # se reads, BC is in column 12
                layout=="PAIRED" & sequencer=="NextSeq",
                "/groups/stark/vloubiere/projects/ORFTRAP_1/git_orftrap_1/functions/demultiplex_pe_12.pl", # pe reads, BC is in column 12 (typically what we use)
                layout=="SINGLE" & sequencer=="HiSeq",
                "/groups/stark/vloubiere/projects/ORFTRAP_1/git_orftrap_1/functions/demultiplex_se_14.pl") # pe reads, BC is in column 12
    cmd <- paste("samtools view -@", cores-1, bam_path, 
                 "| perl", 
                 normalizePath(.f),
                 BC, 
                 fq_prefix, 
                 "; gzip", 
                 paste0(fq_prefix, "_1.fq"))
    if(!is.na(fq2))
      cmd <- paste0(cmd, "; gzip ", fq_prefix, "_2.fq")
    cmd
  }
}, .(layout, barcodes, fq1, fq2)]

## Trim reads
meta[, trim_cmd:= {
  if(!file.exists(fq1_trim))
  {
    cmd <- paste0("trim_galore --gzip -o ", dirname(fq1_trim), "/ ", fq1)
    cmd
  }
}, fq1_trim]

## Alignment ----
meta[, aln_cmd:= {
  if(!file.exists(bam))
  {
    #### BOWTIE 2
    x <- switch(species,
                "mouse" = "/groups/stark/vloubiere/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome",
                "human" = "/groups/stark/vloubiere/genomes/Homo_sapiens/hg38/Bowtie2Index/genome")
    cmd <- paste("bowtie2 -p", cores,
                 "-U", paste0(fq1_trim, collapse= ","),
                 "-x", x,
                 "| samtools sort -@", cores-1, "-o", bam)
  }
}, .(bam, species)]

## Get bam stats ----
meta[, stats_cmd:= {
  if(!file.exists(bam_stats))
    paste("samtools stats -@", cores-1, bam, "| grep ^SN | cut -f 2- >", bam_stats)
}, .(bam, bam_stats)]

## Collapsed bam files ----
meta[, collapse_cmd:= {
  if(!file.exists(bam_unique))
    paste("samtools sort -n -@", cores-1, bam, 
          "| samtools fixmate -m - - | samtools sort -@", cores-1, 
          "| samtools markdup -r - - | samtools view -q 30 -b -o",  bam_unique)
}, .(bam, bam_unique)]

## Compute insertions ----
meta[, insertions_cmd:= {
  if(any(!file.exists(c(counts_same_strand, counts_rev_strand, bed_file))))
  {
    x <- switch(species,
                "mouse" = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_mm10.gtf",
                "human" = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_hg38.gtf")
    paste("Rscript /groups/stark/vloubiere/projects/ORFTRAP_1/git_orftrap_1/functions/bamToBed_and_assign_insertions.R", 
          bam_unique,
          x,
          bed_file,
          gsub("_same_strand.txt$", "", counts_same_strand))
  }
}, .(species, bam_unique, bed_file, counts_same_strand, counts_rev_strand)]

# Submit ----
cols <- c("demultiplex_cmd", "trim_cmd", "aln_cmd", "stats_cmd", "collapse_cmd", "insertions_cmd")
cols <- cols[cols %in% names(meta)]
if(length(cols))
{
  missing_files <- meta[, apply(.SD, 1, function(x) any(!is.na(x))), (meta), .SDcols= cols]$V1
  if(any(missing_files))
  {
    run <- meta[(missing_files)]
    run[, cmd:= paste(na.omit(unique(unlist(.SD))), collapse= ";"), sampleID, .SDcols= patterns("_cmd$")]
    run[, {
      vl_bsub(cmd, 
              cores= cores, 
              m = mem, 
              name = "vlloub", 
              t = '1-00:00:00',
              o= "/groups/stark/vloubiere/projects/ORFTRAP_1/db/logs/",
              e= "/groups/stark/vloubiere/projects/ORFTRAP_1/db/logs/")
    }, cmd]
  }
}

# Save processed metadata ----
meta <- meta[, grep("cmd$", invert = T, names(meta), value = T), with= F]
saveRDS(meta, 
        "Rdata/processed_metadata_ORFtrap.rds")

# Save selected screens ----
sel <- meta[paper_screen!="hRepressor" & !grepl("frame", screen)]
sel[, c("screen", "condition"):= .(paper_screen, paper_condition)]
sel$paper_screen <- sel$paper_condition <- NULL
sel[, Cc:= switch(screen,
                  "Activator"= "#7EBC87",
                  "Repressor"= "#1D71B8",
                  "PTGR"= "#FDC200",
                  "hActivator"= "springgreen3",
                  "hRepressor"= "mediumpurple3"), screen]
saveRDS(sel,
        "Rdata/selected_screens_processed_metadata_ORFtrap.rds")
