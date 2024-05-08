#' Title
#' 
#' @description Function to call hits from a given screen. Receives sorted and unsorted counts as input, computes FC table 
#' 
#' @param sorted_forward_counts Sorted forward counts file (see bamToBed_and_assign_insertions.R function for further details)
#' @param sorted_reverse_counts Sorted reverse counts file (see bamToBed_and_assign_insertions.R function for further details)
#' @param unsorted_forward_counts Unsorted (input) forward counts file (see bamToBed_and_assign_insertions.R function for further details)
#' @param unsorted_reverse_counts Unsorted (input) reverse counts file (see bamToBed_and_assign_insertions.R function for further details)
#' @param exons_start_gtf gtf exon file. Used for consistent ordering of output FC table
#' @param name Name to be appended at the beginning of output file
#' @param output_suffix Suffix to be appended at the end of output file. Default to "_vs_revStrand.txt"
#' @param output_folder_FC_file Output folder for FC files
#'
#' @return Returns FC tables containing DESeq2-like columns
#' @export
#'
#' @examples
#' ORFtrap_call_hits_strandBias(sorted_forward_counts = c("db/gene_assignment/Activator1_sort_rep1_same_strand.txt",
#'                                                        "db/gene_assignment/Activator1_sort_rep2_same_strand.txt"),
#'                              sorted_reverse_counts = c("db/gene_assignment/Activator1_sort_rep1_rev_strand.txt",
#'                                                        "db/gene_assignment/Activator1_sort_rep2_rev_strand.txt"),
#'                              unsorted_forward_counts = c("db/gene_assignment/Activator1_input_rep1_same_strand.txt",
#'                                                          "db/gene_assignment/Activator1_input_rep2_same_strand.txt",
#'                                                          "db/gene_assignment/Repressor2_input_rep1_same_strand.txt",
#'                                                          "db/gene_assignment/Repressor2_input_rep2_same_strand.txt",
#'                                                          "db/gene_assignment/PTGR_input_rep1_same_strand.txt"),
#'                              unsorted_reverse_counts = c("db/gene_assignment/Activator1_input_rep1_rev_strand.txt",
#'                                                          "db/gene_assignment/Activator1_input_rep2_rev_strand.txt",
#'                                                          "db/gene_assignment/Repressor2_input_rep1_rev_strand.txt",
#'                                                          "db/gene_assignment/Repressor2_input_rep2_rev_strand.txt",
#'                                                          "db/gene_assignment/PTGR_input_rep1_rev_strand.txt"),
#'                              exons_start_gtf = "/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/exons_start_mm10.gtf",
#'                              name = "test",
#'                              output_suffix = "_vs_revStrand.txt",
#'                              output_folder_FC_file = "db/FC_tables/")

ORFtrap_call_hits_strandBias <- function(sorted_forward_counts, 
                                         sorted_reverse_counts, 
                                         unsorted_forward_counts,
                                         unsorted_reverse_counts,
                                         exons_start_gtf,
                                         name,
                                         output_suffix= "_vs_revStrand",
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
  inputFw <- rbindlist(lapply(unsorted_forward_counts, fread))
  inputRev <- rbindlist(lapply(unsorted_reverse_counts, fread))
  sampleFw <- rbindlist(lapply(sorted_forward_counts, fread))
  sampleRev <- rbindlist(lapply(sorted_reverse_counts, fread))
  
  dat <- rbindlist(list(count.sample.fw= sampleFw,
                        count.sample.rev= sampleRev,
                        count.input.fw= inputFw,
                        count.input.rev= inputRev),
                   idcol= "cdition")
  dat <- na.omit(dat[dist<5e4])
  dat <- dat[, .(count= .N), .(gene_id, cdition)]
  dat <- dcast(dat,
               gene_id~cdition,
               value.var = "count",
               fill= 0)
  
  # Import gene exons
  genes <- rtracklayer::import(exons_start_gtf)
  genes <- as.data.table(mcols(genes)[, c("gene_id", "gene_name")])
  genes <- unique(genes)
  setorderv(genes, "gene_id")
  
  # Format dat
  dat <- merge(genes, dat, by= "gene_id", all.x= T)
  
  # Fisher test fw vs rev
  dat[count.sample.fw>=3, c("OR", "pval"):= {
    .t <- matrix(c(count.sample.fw+1,
                   count.sample.rev+1,
                   count.input.fw+1,
                   count.input.rev+1),
                 ncol= 2)
    fisher.test(.t, alternative = "greater")[c("estimate", "p.value")]
  }, .(count.sample.fw, count.sample.rev, count.input.fw, count.input.rev)]
  dat[count.sample.fw>=3, log2OR:= log2(OR)]
  dat[count.sample.fw>=3, padj:= p.adjust(pval, "fdr")]
  dat[, hit:= padj<0.05 & log2OR>=1]
  
  # Add binomial test
  dat[count.sample.fw>=3, binom_pval:= {
    binom.test(x = count.sample.fw+1,
               n = count.sample.fw+count.sample.rev+1,
               p = (count.input.fw+1)/(count.input.fw+count.input.rev+1),
               alternative = "greater")["p.value"]
  }, .(count.sample.fw, count.sample.rev, count.input.fw, count.input.rev)]
  dat[count.sample.fw>=3, binom_padj:= p.adjust(binom_pval, "fdr")]
  dat[, binom_hit:= binom_padj<0.001]
  
  # Clean and save
  dat <- genes[dat, on= c("gene_id", "gene_name")]
  dat$OR <- dat$pval <- dat$binom_pval <- NULL
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
