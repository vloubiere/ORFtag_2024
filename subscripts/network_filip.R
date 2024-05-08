setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(data.table)
require(vlfunctions)
require(parallel)
require(GO.db)

#------------------------------------#
# Format data and retrieve complexes
#------------------------------------#
dat <- readRDS("/groups/stark/nemcko/work/0method.paper/analysis/annotation_171hits_all.RDS")
dat <- as.data.table(dat)

#------------------------------------#
# Network
#------------------------------------#
set.seed(5)
net <- vl_STRING_interaction(symbols = dat$hits_all,
                             species = "Mm",
                             plot = F,
                             version= "11.0")

.g1 <- plot(net, score_cutoff = 900)
louvain1 <- igraph::cluster_louvain(.g1)
Cc1 <- rainbow(max(louvain1$membership))

.g2 <- plot(net, score_cutoff = 950)
louvain2 <- igraph::cluster_louvain(.g2)
Cc2 <- rainbow(max(louvain2$membership))

pdf("pdf/hits_string_network_filip.pdf", 10, 10)
# extract igraph
plot(.g1, 
     vertex.color= adjustcolor(Cc1[louvain1$membership], 0.6))
plot(.g2, 
     vertex.color= adjustcolor(Cc2[louvain2$membership], 0.6))
dev.off()

#------------------------------------#
# GO enrichment
#------------------------------------#
GO <- select(org.Mm.eg.db::org.Mm.eg.db,
             keys = unique(dat$geneID_corrected),
             columns = "GO",
             keytype = "ENSEMBL")
GO <- as.data.table(GO)
GO <- na.omit(GO[ENSEMBL %in% dat$geneID_corrected])
GO[data.table(ENSEMBL= dat[louvain1$names, geneID_corrected, on="hits_all"], cl= louvain1$membership), cl:= i.cl, on= "ENSEMBL"]
GO[, check:= any(!is.na(cl)), GO]
GO <- GO[(check), !"check"]
all_genes <- unique(GO$ENSEMBL)
res <- GO[, {
  GO_genes <- ENSEMBL
  .SD[!is.na(cl), {
    tab <- table(factor(all_genes %chin% GO_genes, levels = c(T,F)),
                 factor(all_genes %chin% ENSEMBL, levels = c(T,F)))
    .f <- fisher.test(tab, alternative = "greater")
    .(OR= .f[["estimate"]],
      pval= .f[["p.value"]])
  }, cl]
}, GO]
res[, padj:= -log10(p.adjust(pval, method= "fdr"))]
res[, log2OR:= log2(OR)]
setorderv(res, c("cl", "padj"), c(1, -1))
top <- res[, head(.SD, n= 5), cl]
descript <- as.list(GOTERM)
top$description <- sapply(top$GO, function(x) descript[[x]]@Term)

pdf("pdf/hits_top_go_louvain.pdf", width = 13, height = 30)
my_table_theme <- gridExtra::ttheme_default(core= list(bg_params = list(fill= adjustcolor(Cc1[top$cl], 0.6), col= NA)))
gridExtra::grid.table(top[, !c("GO", "OR", "pval")], theme= my_table_theme)
dev.off()
