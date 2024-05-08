setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(data.table)

# Import human genes ----
human <- rtracklayer::import("db/gtf/unique_protein_coding_transcripts_hg38.gtf")
human <- as.data.table(human)
human <- unique(human[, .(gene_id)])

# Add entrez id ----
hs <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                            keys = human$gene_id,
                            columns = c("ENTREZID", "SYMBOL"),
                            keytype = "ENSEMBL")
hs <- unique(na.omit(as.data.table(hs)))
setnames(hs, c("gene_id", "entrez_id", "symbol"))
human <- merge(human,
               hs,
               by= "gene_id",
               all.x= T)

# Retrieve mouse/human homologs ----
# download 25/11/21: "http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt")
mgi <- fread("db/public_db/HOM_MouseHumanSequence_20230109.rpt")
hs <- mgi[`Common Organism Name`=="human", 
          .(`DB Class Key`,
            entrez_id= as.character(`EntrezGene ID`))]
human[hs, homolog_id:= `i.DB Class Key`, on= "entrez_id"]
mm <- mgi[`Common Organism Name`=="mouse, laboratory", 
          .(homolog_id= `DB Class Key`,
            mouse_entrez= as.character(`EntrezGene ID`))]
human[mm, mouse_entrez:= `i.mouse_entrez`, on= "homolog_id"]

# Add mouse ENSEMBL IDs ----
extra <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
                               keys = human$mouse_entrez,
                               columns = "ENSEMBL",
                               keytype = "ENTREZID")
extra <- unique(na.omit(as.data.table(extra)))
human[extra, mouse_ENSEMBL:= i.ENSEMBL, on= "mouse_entrez==ENTREZID"]

# fusion proteins, v97 on Jan. 2023 ----
cosmic <- fread("/groups/stark/nemcko/ORFtrap_paper/db/CosmicFusionExport_v97.tsv.gz")
cosmic <- as.data.table(cosmic)
cosmic[, `5'_GENE_NAME`:= gsub("(.+)_.+", "\\1", `5'_GENE_NAME`)]
cosmic[, `3'_GENE_NAME`:= gsub("(.+)_.+", "\\1", `3'_GENE_NAME`)]
human[, fusion:= symbol %in% cosmic[, c(`5'_GENE_NAME`, `3'_GENE_NAME`)]]

# TF genes ----
tf <- readxl::read_xlsx("/groups/stark/nemcko/work/ORFtrap/general/Lambert-el-al.(2018).xlsx",
                        sheet = 2,
                        skip = 1)
tf <- as.data.table(tf[, c(1,3,4)])
setnames(tf, c("gene_id", "DBD", "TF"))
human[, TF:= gene_id %in% tf[TF=="Yes", gene_id]]

# Activation/Repression Domains ----
soto <- readxl::read_xlsx("db/public_db/Soto_2022.xlsx", sheet = 5)
soto <- as.data.table(soto)
soto <- soto[`N or S?` %in% c("S", "N and S") # Necessary or Necessary and Sufficient
             & `Confidence (H, M or L)`=="H"] # High confidence
soto[, domain_type:= fcase(all(`Domain type`=="AD"), "AD", # For each gene, check if it contains Activating Domains (AD), Repressive Domains (RD) or both (Bif for Bifunctional) 
                           all(`Domain type`=="RD"), "RD",
                           default = "Bif"), `ENSEMBL gene ID`]
human[, AD:= gene_id %in% soto[domain_type=="AD", `ENSEMBL gene ID`]]
human[, RD:= gene_id %in% soto[domain_type=="RD", `ENSEMBL gene ID`]]

# RNA binding Domains ----
# "https://rbpbase.shiny.embl.de/data/RBPbase_Hs_DescriptiveID.xlsx"
RBP <- readxl::read_xlsx("db/public_db/RBPbase_Hs_DescriptiveID.xlsx", col_types = "text")
RBP <- as.data.table(RBP)
RBP <- RBP[any_Hs=="YES" & `Gene-type-Hs\nRBPANNO000000065.1`=="protein_coding" & hits_Hs>1]
RBP <- unique(unlist(strsplit(RBP$ID, "|", fixed = T)))
human[, RBP:= gene_id %in% RBP]

# Alerasool hits ----
al <- readxl::read_xlsx("db/public_db/alerasool_2022.xlsx", sheet = 2)
al <- as.data.table(al)
human[, Alerasool:= entrez_id %in% al[Hit=="hit", `Gene ID`]]

# DelRosso hits ----
DelRosso <- readxl::read_xlsx("db/public_db/DelRosso_sup_Table_1.xlsx", sheet = 1)
DelRosso <- as.data.table(DelRosso)
human[, DelRosso:= gene_id %in% DelRosso[`Hit minCMV`=="Hit" & !is.na(`Hit minCMV`), `Ensembl ID`]]

# Save ----
saveRDS(human,
        "Rdata/Annotations_hs_homologs.rds")
