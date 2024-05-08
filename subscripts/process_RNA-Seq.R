setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(data.table)
require(parallel)
require(vlfunctions)

# Import metadata ----
meta <- readRDS("Rdata/RNA_trimmed_metadata.rds")

# Parse Metadata ----
if(!all(meta[, sampleID==paste0(screen, "_", condition, "_rep", replicate)]))
  stop("sampleID should be the catenation of screen, condtion and replicate, so that re-sequenced samples are correctly collapsed")

# Generate output paths ----
meta[, bam:= paste0("/groups/stark/vloubiere/projects/ORFTRAP_1/db/bam/", sampleID, "_RNA.bam")]
meta[, bed:= paste0("/groups/stark/vloubiere/projects/ORFTRAP_1/db/bed/", sampleID, "_RNA.bed")]
meta[, exon_assignment:= paste0("/groups/stark/vloubiere/projects/ORFTRAP_1/db/exon_assignment/", sampleID, "_RNA.rds")]
meta[, bw:= paste0("/groups/stark/vloubiere/projects/ORFTRAP_1/db/bw/", sampleID, "_RNA.bw")]

# Save processed metadata ----
saveRDS(meta,
        "Rdata/RNA_processed_metadata.rds")

# mm10 index
if(!file.exists("/groups/stark/vloubiere/genomes/Mus_musculus/subreadr_mm10/subreadr_mm10_index.log"))
  Rsubread::buildindex(basename= "/groups/stark/vloubiere/genomes/Mus_musculus/subreadr_mm10/subreadr_mm10_index",
                       reference= "/groups/stark/vloubiere/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa")

# Alignment ----
meta[, {
  if(!file.exists(bam))
  {
    print(paste0("START...", bam))
    Rsubread::align(index= "/groups/stark/vloubiere/genomes/Mus_musculus/subreadr_mm10/subreadr_mm10_index",
                    readfile1= fq1_frame,
                    # readfile2= fq2,
                    input_format = "gzFASTQ",
                    maxMismatches = 2,
                    nthreads = 8,
                    unique = T,
                    output_file= bam)
  }
  print(paste(bam, "--> DONE!"))
}, .(fq1_frame, bam)]

# Load the modules that will be used by the pipeline ----
meta[, load_cmd:= paste0(c("module load build-env/2020",
                           "module load trim_galore/0.6.0-foss-2018b-python-2.7.15",
                           "module load cutadapt/1.18-foss-2018b-python-2.7.15",
                           "module load bowtie2/2.3.4.2-foss-2018b",
                           "module load samtools/1.9-foss-2018b",
                           "module load r/3.6.2-foss-2018b"), collapse = "; ")]

# Collapse unique reads and extract frame ----
meta[, uniq_bed_cmd:= {
  if(!file.exists(bed))
    paste("Rscript /groups/stark/vloubiere/projects/ORFTRAP_1/git_orftrap_1/functions/collapse_RNASeq_reads.R", 
          bam,
          bed)
}, .(bam, bed)]

# Assign reads to genomic features ----
meta[, assign_cmd:= {
  if(!file.exists(exon_assignment))
    paste("Rscript /groups/stark/vloubiere/projects/ORFTRAP_1/git_orftrap_1/functions/compute_frame_RNA_reads.R", 
          bed,
          "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_phase_mm10.gtf",
          "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/transcripts_and_upstream_region_mm10.gtf",
          exon_assignment)
}, .(bed, exon_assignment)]

# Submit jobs ----
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
            cores= 6, 
            m = 20, 
            name = "vlloub", 
            t = '04:00:00',
            o= "/groups/stark/vloubiere/projects/ORFTRAP_1/db/logs/",
            e= "/groups/stark/vloubiere/projects/ORFTRAP_1/db/logs/")
  }, .(cmd= V1)]
}

stop()

# Generate bw file ----
meta[, {
  if(!file.exists(bw))
  {
    .c <- rtracklayer::import.bed(bed)
    cov <- coverage(.c)/length(.c)*1e6
    rtracklayer::export.bw(GRanges(cov),
                           con= bw)
  }
  print("DONE")
}, .(bed, bw)]