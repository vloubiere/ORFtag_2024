setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
source("git_orftrap_1/functions/ORFTRAP_call_hits.R")
source("git_orftrap_1/functions/ORFTRAP_call_hits_strandBias.R")

meta <- readRDS("Rdata/selected_screens_processed_metadata_ORFtrap.rds")
meta <- unique(meta[, .(species, screen, condition, counts_same_strand, counts_rev_strand)])
mouse <- meta[species=="mouse"]
human <- meta[species=="human"]
frame <- readRDS("Rdata/processed_metadata_ORFtrap.rds")[grepl("frame", screen)]

# Mouse screens: using merged Input ----
mouse[condition=="sort", {
  ORFtrap_call_hits(sorted_forward_counts = counts_same_strand,
                    unsorted_forward_counts = mouse[condition=="input", counts_same_strand],
                    exons_start_gtf = "db/gtf/exons_start_mm10.gtf", 
                    name = screen,
                    output_suffix = "_vs_allMergedInputs.txt", 
                    output_folder_FC_file = normalizePath("db/FC_tables/selected_screens", mustWork = F))
}, screen]

# Mouse screens: using reversed insertions from merged Input ----
mouse[condition=="sort", {
  ORFtrap_call_hits(sorted_forward_counts = counts_rev_strand,
                    unsorted_forward_counts = mouse[condition=="input", counts_rev_strand],
                    exons_start_gtf = "db/gtf/exons_start_mm10.gtf", 
                    name = screen,
                    output_suffix = "_rev_counts_vs_allMergedInputs.txt", 
                    output_folder_FC_file = normalizePath("db/FC_tables/selected_screens", mustWork = F))
}, screen]

# Human Activator screen: using merged Input ----
human[condition=="sort", {
  ORFtrap_call_hits(sorted_forward_counts = counts_same_strand, 
                    unsorted_forward_counts = human[condition=="input", counts_same_strand],
                    exons_start_gtf = "db/gtf/exons_start_hg38.gtf",
                    name = screen,
                    output_suffix = "_vs_allMergedInputs.txt",
                    output_folder_FC_file = normalizePath("db/FC_tables/selected_screens", mustWork = F))
}, screen]

# Human Activator screen: using reversed insertions from merged Input ----
human[condition=="sort", {
  ORFtrap_call_hits(sorted_forward_counts = counts_rev_strand, 
                    unsorted_forward_counts = human[condition=="input", counts_rev_strand],
                    exons_start_gtf = "db/gtf/exons_start_hg38.gtf",
                    name = screen,
                    output_suffix = "_rev_counts_vs_allMergedInputs.txt",
                    output_folder_FC_file = normalizePath("db/FC_tables/selected_screens", mustWork = F))
}, screen]


# Frame-specific screens: using individual inputs ----
frame[, {
  ORFtrap_call_hits(sorted_forward_counts= counts_same_strand[condition=="sort"],
                    unsorted_forward_counts= counts_same_strand[condition=="input"],
                    exons_start_gtf = "db/gtf/exons_start_mm10.gtf", 
                    name = screen,
                    output_suffix = "_vs_frameSpecificInput.txt", 
                    output_folder_FC_file = normalizePath("db/FC_tables/selected_screens", mustWork = F))
}, screen]

# Frame-specific screens: using merged input from previous screens (for comparison) ----
ORFtrap_call_hits(sorted_forward_counts= frame[condition=="sort", counts_same_strand],
                  unsorted_forward_counts= mouse[condition=="input", counts_same_strand],
                  exons_start_gtf = "db/gtf/exons_start_mm10.gtf", 
                  name = "three_frames_merged",
                  output_suffix = "_vs_allMergedInputs.txt", 
                  output_folder_FC_file = normalizePath("db/FC_tables/selected_screens", mustWork = F))
