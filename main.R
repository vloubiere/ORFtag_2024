setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")

if(F)
{
  # Custom annotation files ----
  file.edit("ORFtag_2024/subscripts/compute_gtf_non_first_exons_gtf.R")
  file.edit("ORFtag_2024/subscripts/compute_gtf_exons_with_CDS_frame.R")
  file.edit("ORFtag_2024/subscripts/compute_gtf_transcripts_and_upstream_regions.R")
  file.edit("ORFtag_2024/subscripts/compute_gtf_unique_protein_coding_transcripts.R")
  # Parse public data
  file.edit("ORFtag_2024/subscripts/split_full_uniprot_table.R") # Get protein domains info
  # Human homologs data
  file.edit("ORFtag_2024/subscripts/retrieve_human_databases.R") # TF, RNA-binding proteins...
  
  # Functions ----
  # Demultiplex
  file.edit("ORFtag_2024/functions/demultiplex_pe_12.pl")
  file.edit("ORFtag_2024/functions/demultiplex_se_12.pl")
  file.edit("ORFtag_2024/functions/demultiplex_pe_14.pl")
  file.edit("ORFtag_2024/functions/demultiplex_se_14.pl")
  # Assign insertions to exons
  file.edit("ORFtag_2024/functions/bamToBed_and_assign_insertions.R")
  # Calling ORFtag hits
  file.edit("ORFtag_2024/functions/ORFTRAP_call_hits.R")
  file.edit("ORFtag_2024/functions/ORFTRAP_call_hits_strandBias.R")
  # GO terms
  file.edit("ORFtag_2024/functions/GO_enrichments.R")
  # RNA-Seq
  file.edit("ORFtag_2024/functions/process_ORFtag_RNASeq_reads.pl")
  file.edit("ORFtag_2024/functions/collapse_RNASeq_reads.R")
  file.edit("ORFtag_2024/functions/compute_frame_RNA_reads.R")
  
  # Data processing ----
  # ORFtag pipeline
  file.edit("ORFtag_2024/subscripts/ORFTRAP_pipeline.R")
  file.edit("ORFtag_2024/subscripts/generate_bw_files.R")
  # Calling ORFtag hits
  file.edit("ORFtag_2024/subscripts/alignment_statistics.R")
  file.edit("ORFtag_2024/subscripts/ORFTRAP_enrichment.R")
  # Final tables
  file.edit("ORFtag_2024/subscripts/Final_tables.R")
  # RNA-Seq screens (transcripts frame)
  file.edit("ORFtag_2024/subscripts/trim_VBC_reads.R")
  file.edit("ORFtag_2024/subscripts/process_RNA-Seq.R")
  
  # Figures ----
  # Figure 1
  file.edit("ORFtag_2024/subscripts/volcanos_screens.R") # Hits volcano plot for each screen
  # Figure 2
  file.edit("ORFtag_2024/subscripts/overlaps_screens.R") # Veen diag overlap between the 3 screens
  file.edit("ORFtag_2024/subscripts/human_annotations_enrichments.R")# Human TF genes and domains 
  file.edit("ORFtag_2024/subscripts/GO_enrichment_mouse_hits.R")# GO enrich
  file.edit("ORFtag_2024/subscripts/boxplots_validations.R") # Candidates validations FACS
  # Figure 3
  file.edit("ORFtag_2024/subscripts/average_tracks_input.R") # Average tracks
  file.edit("ORFtag_2024/subscripts/saturation_input_screens.R") # Fraction of genes being trapped
  file.edit("ORFtag_2024/subscripts/barplot_trapped_genes_hits_CDS_length.R") # Ratio trapped genes CDS including hits
  file.edit("ORFtag_2024/subscripts/barplot_trapped_genes_TPM_quantiles.R") # % trapped genes CDS and TPM quantiles
  
  # Supplementary Figures ----
  # Supplementary Figure 2
  file.edit("ORFtag_2024/subscripts/review_compute_insertions_exons_assignment.R") # Frame-specific screens
  file.edit("ORFtag_2024/subscripts/review_compare_frame_specific_screens_mouse.R")
  file.edit("ORFtag_2024/subscripts/FACS_Dox_induced.R") # FACS DOX induced TetOFF
  file.edit("ORFtag_2024/subscripts/review_RNASeq_reads_frame_refined.R")
  # Supplementary Figure 3
  file.edit("ORFtag_2024/subscripts/mouse_STRING_networks.R") # String protein interaction networks hits
  file.edit("ORFtag_2024/subscripts/PCC_replicates.R") # Correlations between replicates
  file.edit("ORFtag_2024/subscripts/scatterplots_comparisons_screens.R") # Compare FDRs between screens
  file.edit("ORFtag_2024/subscripts/review_enrichments_Alerasool_DelRosso.R") # Comparison hits/Alerasool enrich for TF/AD proteins
  # Supplementary Figure 4
  file.edit("ORFtag_2024/subscripts/review_first_CDS_analysis.R") # Missing fraction of tagged prot
  file.edit("ORFtag_2024/subscripts/first_exons_analysis.R") # Check if first exon contains coding seq
  file.edit("ORFtag_2024/subscripts/intronless_genes_analysis.R") # Intronless genes protein domain enrichment
}