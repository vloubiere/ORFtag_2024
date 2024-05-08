#' Title
#' 
#' @description Function to call hits from a given screen. Receives sorted and unsorted counts as input, computes FC table 
#' 
#' @param sorted_forward_counts Sorted forward counts file (see bamToBed_and_assign_insertions.R function for further details)
#' @param unsorted_forward_counts Unsorted (input) forward counts file (see bamToBed_and_assign_insertions.R function for further details)
#' @param exons_start_gtf gtf exon file. Used for consistent ordering of output FC table
#' @param name Name to be appended at the beginning of output file
#' @param output_suffix Suffix to be appended at the end of output file. Default to "_vs_unsorted.txt"
#' @param output_folder_FC_file Output folder for FC files
#'
#' @return Returns FC tables containing DESeq2-like columns
#' @export
#'
#' @examples
#' ORFtrap_call_hits(sorted_forward_counts = c("db/gene_assignment/Activator1_sort_rep1_same_strand.txt",
#'                                             "db/gene_assignment/Activator1_sort_rep2_same_strand.txt"),
#' unsorted_forward_counts = c("db/gene_assignment/Activator1_input_rep1_same_strand.txt",
#'                             "db/gene_assignment/Activator1_input_rep2_same_strand.txt",
#'                             "db/gene_assignment/Repressor2_input_rep1_same_strand.txt",
#'                             "db/gene_assignment/Repressor2_input_rep2_same_strand.txt",
#'                             "db/gene_assignment/PTGR_input_rep1_same_strand.txt"),
#' exons_start_gtf = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_mm10.gtf",
#' name = "test",
#' output_suffix = "_vs_unsorted.txt",
#' output_folder_FC_file = "db/FC_tables/")

ORFtrap_call_hits <- function(sorted_forward_counts, 
                              unsorted_forward_counts,
                              exons_start_gtf,
                              name,
                              output_suffix= "_vs_unsorted.txt",
                              output_folder_FC_file)
{
  require(rtracklayer)
  require(data.table)
  require(GenomicRanges)
  
  # Checks
  if(!is.character(name) | length(name)!=1)
    stop("name should be a unique character specifying the name of the screen")
  if(anyDuplicated(sorted_forward_counts))
    stop("Duplicated filenames in sorted_forward_counts")
  if(anyDuplicated(unsorted_forward_counts))
    stop("Duplicated filenames in unsorted_forward_counts")
  if(!is.character(output_suffix) | length(output_suffix)!=1)
    stop("output_suffix should be a unique character specifying the name of the screen")
  
  # Import counts
  input <- rbindlist(lapply(unsorted_forward_counts, fread))
  sample <- rbindlist(lapply(sorted_forward_counts, fread))
  dat <- rbindlist(list(count.input= input, count.sample= sample), idcol= "cdition")
  dat <- na.omit(dat[dist<2e5])
  total <- dat[, .N, .(cdition)]
  dat <- dat[, .(count= .N), .(gene_id, cdition)]
  dat <- dcast(dat, gene_id~cdition, value.var = "count")
  
  # Import gene exons
  genes <- rtracklayer::import(exons_start_gtf)
  genes <- as.data.table(mcols(genes)[, c("gene_id", "gene_name")])
  genes <- unique(genes)
  setorderv(genes, "gene_id")
  
  # Format dat
  dat <- merge(genes, dat, by= "gene_id", all.x= T)
  dat[is.na(count.input), count.input:= 0]
  dat[is.na(count.sample), count.sample:= 0]
  dat[, total.input:= total[cdition=="count.input", N]]
  dat[, total.sample:= total[cdition=="count.sample", N]]
  
  # Fisher test sample vs input
  dat[count.sample>=3, c("OR", "pval"):= {
    .t <- matrix(c(total.input-count.input+1,
                   count.input+1,
                   total.sample-count.sample+1,
                   count.sample+1),
                 ncol= 2)
    fisher.test(.t, alternative = "greater")[c("estimate", "p.value")]
  }, .(count.input, total.input, count.sample, total.sample)]
  dat[count.sample>=3, log2OR:= log2(OR)]
  dat[count.sample>=3, padj:= p.adjust(pval, "fdr")]
  dat[, hit:= padj<0.001 & log2OR>=1]
  
  # Clean and save
  dat <- genes[dat, on= c("gene_id", "gene_name")]
  dat$OR <- dat$pval <- NULL
  dir.create(output_folder_FC_file, showWarnings = F, recursive = T)
  FC_table <- paste0(output_folder_FC_file, "/", name, output_suffix)
  fwrite(dat, 
         FC_table, 
         col.names = T, 
         row.names = F, 
         sep= "\t",
         quote= F, 
         na= NA)
  return(paste0(name, ": ", sum(dat$hit, na.rm = T), " hits were called!\nFC file -> ", FC_table, "\n"))
}
