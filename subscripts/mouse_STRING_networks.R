setwd("/groups/stark/vloubiere/projects/ORFTRAP_1/")
require(vlfunctions)
require(igraph)
require(ggplot2)
require(ggraph)

# Import data ----
annot <- readRDS("Rdata/Annotations_hs_homologs.rds")
mouse <- readRDS("Rdata/final_table_mouse.rds")
human <- readRDS("Rdata/final_table_human.rds")

# Extract human hits ----
human_hits <- human[(hit)]
human_hits <- annot[human_hits, .(screen, mouse_ENSEMBL), on= "human_ENSEMBL==gene_id"]
human_hits <- na.omit(human_hits)

# Compute classes in mouse ----
mouse[human_hits, Cc:= fcase(gene_id %in% i.mouse_ENSEMBL, "tomato", # Found in human
                             gene_id %in% annot$mouse_ENSEMBL, "cornflowerblue", # Homolog but not found
                             default = "grey"), on= "screen"] # No Homolog

# Import mouse database and compute networks ----
STRINGdb <- vl_STRING_getDB(species = "Mm",
                            network_type = "full",
                            version= "11.0")

# Random sampling ----
pdf("pdf/mouse_hits_string_ppi_likelyhood.pdf", 5, 5)
par(mfrow= c(2,2),
    tcl= -0.2,
    mgp= c(2,0.5, 0),
    las= 1)
mouse[, number_hits:= sum(hit, na.rm= T), screen]
mouse[, {
  # Extract all interactions > 900
  obj <- vl_STRING_interaction(STRINGdb = STRINGdb,
                               symbols = gene_name,
                               score.cutoff = 900,
                               remove.non.connected = T,
                               plot= F)
  edges <- obj$E
  # Count the number of hits with at least one other interaction with another hit 
  hits.i <- length(unique(edges[from %in% gene_name[(hit)] & to %in% gene_name[(hit)], c(from, to)]))
  # Randomly sample 100 times and count number of proteins with at least 1 interaction
  rdm.i <- sapply(1:100, function(x)
  {
    set.seed(x)
    .s <- sample(gene_name, number_hits)
    length(unique(edges[from %in% .s & to %in% .s, c(from, to)]))
  })
  # Histogram
  rdm_frac <- rdm.i/number_hits
  hits_frac <- hits.i/number_hits
  hist(rdm_frac,
       freq = F,
       xlab= "Fraction of interacting proteins",
       xlim= c(0, .5),
       main= screen)
  abline(v= hits_frac, lwd= 2)
  print(screen)
}, .(screen, number_hits)]
dev.off()

# Plot networks ----
mouse[(hit), size:= (log2OR-min(log2OR, na.rm = T))/diff(range(log2OR, na.rm = T))*8+0.5, screen]

# Plot and color based on louvain ----
pdf("pdf/mouse_hits_string_network_per_screen_louvain_coloring.pdf", 3.75, 3.75)
par(las= 1,
    tcl= -0.2,
    mgp= c(2, 0.5, 0),
    mar= c(8,8,8,8))
mouse[(hit), {
  # Extract interactions
  obj <- vl_STRING_interaction(STRINGdb = STRINGdb,
                               symbols = gene_name,
                               score.cutoff = 900,
                               remove.non.connected = T,
                               plot= F,
                               size= size,
                               col= Cc,
                               label.cex= size/10,
                               frame.width = 0.5)
  .i <- vl_STRING_to_igraph(obj)
  # Louvain clustering
  louvain <- igraph::cluster_louvain(.i)
  Cc <- rainbow(max(louvain$membership))
  Cc <- sapply(Cc, function(x) colorRampPalette(c("gray90", x))(3)[2])
  # Color based on Louvain
  set.seed(1)
  pl <- ggraph(.i, layout = "nicely") +
    ggtitle(screen) +
    geom_edge_link0(edge_colour = "grey66") +
    geom_node_point(aes(size = size),
                    fill= Cc[louvain$membership],
                    stroke= 0.1,
                    shape = 21) +
    geom_node_text(aes(label = name,
                       size = size/10),
                   family = "Helvetica",
                   segment.size= 0.2,
                   repel = TRUE,
                   max.overlaps= 10) +
    theme_graph(border= FALSE, base_family = 'Helvetica')+
    theme(legend.position = "none")
  plot(pl)
  print("")
}, screen]
dev.off()

# Plot and color based Ovelrap with Human ----
pdf("pdf/mouse_hits_string_network_per_screen_humanOverlap_coloring.pdf", 3.75, 3.75)
par(las= 1,
    tcl= -0.2,
    mgp= c(2, 0.5, 0),
    mar= c(8,8,8,8))
mouse[(hit) & screen!="PTGR", {
  # Extract interactions
  obj <- vl_STRING_interaction(STRINGdb = STRINGdb,
                               symbols = gene_name,
                               score.cutoff = 900,
                               remove.non.connected = T,
                               plot= F,
                               size= size,
                               col= Cc,
                               label.cex= size/10,
                               frame.width = 0.5)
  .i <- vl_STRING_to_igraph(obj)
  # Color based on overlap with Human
  set.seed(1)
  pl <- ggraph(.i, layout = "nicely") +
    ggtitle(screen) +
    geom_edge_link0(edge_colour = "grey66") +
    geom_node_point(aes(size = size),
                    fill= V(.i)$color,
                    stroke= 0.1,
                    shape = 21) +
    geom_node_text(aes(label = name,
                       size = size/10),
                   family = "Helvetica",
                   segment.size= 0.2,
                   repel = TRUE,
                   max.overlaps= 10) +
    theme_graph(border= FALSE, base_family = 'Helvetica')+
    theme(legend.position = "none")
  plot(pl)
  print("")
}, screen]
dev.off()