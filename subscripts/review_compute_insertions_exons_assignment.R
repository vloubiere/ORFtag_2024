setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(vlfunctions)

meta <- readRDS("Rdata/processed_metadata_ORFtrap.rds")
meta <- meta[grepl("frame", sampleID)]

# Import closest downstream exon assignment ----
dat <- meta[, fread(counts_same_strand), .(screen, condition)]
dat[, screen:= gsub("^Activator2_f", "F", screen)]
dat[, expectedPhase:= switch(screen,
                             "Frame0"= 0L,
                             "Frame1"= 2L,
                             "Frame2"= 1L), screen]

# Add whether exon is in a hit or not
actFrame <- list(Frame0= "db/FC_tables/selected_screens/Activator2_frame0_vs_frameSpecificInput.txt",
                 Frame1= "db/FC_tables/selected_screens/Activator2_frame1_vs_frameSpecificInput.txt",
                 Frame2= "db/FC_tables/selected_screens/Activator2_frame2_vs_frameSpecificInput.txt")
actFrame <- lapply(actFrame, fread)
actFrame <- rbindlist(actFrame,
                      idcol = "screen")
actFrame <- actFrame[(hit)]
dat$hit <- actFrame[dat, .N>0, .EACHI, on= c("screen", "gene_id")]$V1

# Retrieve exon phase
frame <- rtracklayer::import("db/gtf/exons_phase_mm10.gtf")
frame <- as.data.table(frame)
frame[, exonPhase:= as.integer(exonPhase)]
dat$compatiblePhase <- frame[dat, any(exonPhase==expectedPhase), .EACHI, on= "exon_id"]$V1

saveRDS(dat, 
        "db/exon_assignment/Activator2_frame_specific_insertions.rds")
