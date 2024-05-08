GO_enrichment <- function(list_gene_id,
                          species,
                          exon_file,
                          output_file)
{
  # Import
  genes <- rtracklayer::import(exon_file)
  genes <- unique(as.data.table(genes)[, .(gene_id, gene_name)])
  
  # Retrieve annotations ----
  ## Protein domains ----
  domains <- AnnotationDbi::select(switch(species,
                                          "mouse"= EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                                          "human"= EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79),
                                   keys = genes$gene_id, 
                                   keytype = "GENEID", 
                                   columns = "PROTEINDOMAINID")
  domains <- unique(na.omit(as.data.table(domains)))
  setnames(domains, c("gene_id",  "ID"))
  extra <- fread("db/public_db/protein_domains_description_uniqID.txt")
  domains <- merge(domains, 
                   extra, 
                   by.x= "ID", 
                   by.y= "ID")
  ## GO ----
  GO <- AnnotationDbi::select(switch(species,
                                     "mouse"= org.Mm.eg.db::org.Mm.eg.db,
                                     "human"= org.Hs.eg.db::org.Hs.eg.db), 
                              keys = genes$gene_id, 
                              keytype = "ENSEMBL", 
                              columns = "GO")
  GO <- unique(na.omit(as.data.table(GO)))
  GO <- GO[, !"EVIDENCE"]
  setnames(GO, c("gene_id",  "ID", "type"))
  extra <- select(GO.db::GO.db, 
                  keys = GO$ID, 
                  keytype = "GOID", 
                  columns = "TERM")
  extra <- unique(na.omit(as.data.table(extra)))
  setnames(extra, "TERM", "description")
  GO <- merge(GO, 
              extra, 
              by.x= "ID", 
              by.y= "GOID")
  # CC
  CC <- fread("db/public_db/subcellular_location.txt")
  CC <- merge(genes[, .(gene_id)], 
              CC,
              by= "gene_id")
  setnames(CC, "Subcellular_location", "description")
  CC[, type:= "Subcellular_localization"]
  CC[, ID:= NA]
  # Bind all
  annot <- rbindlist(list(Protein_domains= domains,
                          GO= GO,
                          Localization= CC),
                     idcol = "group", 
                     fill= T)
  annot <- unique(annot[, .(group, gene_id, description, type, ID)])
  
  # Enrichment per screen ----
  dat <- lapply(list_gene_id, function(x) copy(genes)[, hit:= gene_id %in% x])
  dat <- rbindlist(dat, idcol = "screen")
  dat[, gene_total:= .N, screen]
  dat[, hit_total:= sum(hit, na.rm = T), screen]
  dat <- merge(dat,
               annot,
               by= "gene_id", allow.cartesian= T)
  dat <- dat[, .(hit_counts= sum(hit, na.rm = T), 
                 annot_counts= .N, 
                 genes= paste0(unique(gene_name[(hit)]), collapse = ",")), 
             .(group, type, screen, description, ID, gene_total, hit_total)]
  dat[, c("estimate", "pval"):= {
    .t <- matrix(c(hit_counts, 
                   hit_total-hit_counts,
                   annot_counts-hit_counts,
                   gene_total-(hit_total+annot_counts-hit_counts)), 
                 ncol= 2)
    fisher.test(.t, alternative = "greater")[c("estimate", "p.value")]
  }, .(gene_total, hit_total, hit_counts, annot_counts)]
  dat[, padj:= p.adjust(pval, "fdr"), .(screen, group, type)]
  dat[, log2OR:= log2(estimate)]
  
  # Select top padj for each description and cap log2OR ----
  dat <- dat[, .SD[which.min(padj)], .(group, type, screen, description)]
  dat[, log2OR:= ifelse(log2OR==Inf, max(log2OR[is.finite(log2OR)]), log2OR), .(group, type, screen)]
  
  # Save ----
  setorderv(dat, c("group", "type", "screen", "log2OR", "padj"), c(1,1,1,-1,1))
  dat$pval <- dat$estimate <- NULL
  dat <- dat[padj<0.05 & hit_counts>=3 & annot_counts>5]
  fwrite(dat,
         output_file,
         col.names = T,
         row.names = F,
         sep= "\t",
         quote= F,
         na= NA)
}