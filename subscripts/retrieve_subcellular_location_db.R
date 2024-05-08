setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(data.table)
require(org.Mm.eg.db)
require(AnnotationDbi)
require(UniProt.ws)
require(tidyverse)
require(biomaRt)

gtf <- rtracklayer::import("db/gtf/exons_start_mm10.gtf")
dat <- as.data.table(gtf)
dat <- dat[, .(gene_id= gsub("[.*]..*", "\\1", gene_id),
               gene_name)]
dat <- unique(dat)

# Add UNIPROT and ENTREZ IDs
extra <- biomaRt::select(org.Mm.eg.db,
                keys = dat$gene_id,
                columns = c("ENTREZID", "UNIPROT"),
                keytype = "ENSEMBL")
extra <- as.data.table(extra)
extra <- unique(na.omit(extra))
setnames(extra, c("gene_id", "entrez_id", "UNIPROT"))
dat <- merge(dat,
             extra,
             by= "gene_id")

# load mouse uniprot dataset from uniprot.ws package
mouseUp <- UniProt.ws(10090)

# load subcellular location for all uniprot IDs
dt <- UniProt.ws::select(mouseUp, 
                         columns = c("cc_subcellular_location"), 
                         keys = dat$UNIPROT, 
                         keytype = "UniProtKB")
dt <- as.data.table(dt)

# Simplify locations
dt[, sub:= tstrsplit(Subcellular.location..CC., " Note=", keep= 1)]
dt[, sub:= gsub("\\s*(\\{[^{}]*(?:(?1)[^{}]*)*\\})", "", sub, perl=TRUE)]
dt[, sub:= tstrsplit(sub, "^SUBCELLULAR LOCATION: ", keep= 2)]
dt[, sub:= gsub("; SUBCELLULAR LOCATION:", "", sub)]
# remove isoforms
dt[, sub:= gsub("\\[", "{", sub)]
dt[, sub:= gsub("\\]:", "}", sub)]
dt[, sub:= gsub("\\s*(\\{[^{}]*(?:(?1)[^{}]*)*\\})", "", sub, perl=TRUE)]
dt[,sub:= gsub("\\.$", "", sub)]
dt[, sub := str_squish(sub)]

final <- dt[, list(sub = unlist(strsplit(sub, "\\."))), by= From]
final <- unique(final)
# remove unnecessary gaps
final[, sub := str_squish(sub)]

#662 unique locations, some too many details, removed everything behind the first ";"
final[, sub2:= tstrsplit(sub, ";", keep= 1)]
final <- final[!is.na(sub2), .(UNIPROT= From, Subcellular_location= sub2)]
final <- merge(dat, 
               final,
               by= "UNIPROT")
final[dat, gene_id:= i.gene_id, on= "UNIPROT"]

# not only 237 locations
fwrite(final,
       "db/public_db/subcellular_location.txt",
       quote= F,
       sep= "\t",
       col.names = T, 
       row.names = F,
       na= NA)