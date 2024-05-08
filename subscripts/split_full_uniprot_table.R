require(data.table)

# Full table was downloaded from "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz"
# Need ~120 G of ram to process it -> downloaded on the scratch
dat <- fread("/scratch/stark/vloubiere/protein2ipr_v92.0.dat.gz")
uniq_IDs <- unique(dat[, .(description= V3, ID= V4, INTERPRO= V2)])

# Add details
download.file("https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list", 
              destfile = "db/public_db/interpro_entry.list")
extra <- fread("db/public_db/interpro_entry.list")
uniq_IDs[extra, type:= i.ENTRY_TYPE, on="INTERPRO==ENTRY_AC"]
fwrite(uniq_IDs,
       "/groups/stark/vloubiere/projects/ORFTRAP_1/db/public_db/protein_domains_description_uniqID.txt",
       col.names = T,
       row.names = F, 
       sep= "\t",
       na= NA,
       quote= F)
