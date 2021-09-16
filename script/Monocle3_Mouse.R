library(monocle3)

data <- as(as.matrix(brain_3prime_bcell[['RNA']]@data), 'sparseMatrix')
pd <- brain_3prime_bcell@meta.data
fd <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
cds <- new_cell_data_set(data, cell_metadata = pd, gene_metadata = fd)

cds@reducedDims$PCA <- brain_3prime_bcell@reductions$pca@cell.embeddings
cds@reducedDims$UMAP <- brain_3prime_bcell@reductions$umap@cell.embeddings

recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

cds@clusters$UMAP$clusters <- brain_3prime_bcell$seurat_clusters
cds@preprocess_aux$gene_loadings <- brain_3prime_bcell@reductions$pca@feature.loadings

cds <- learn_graph(cds)
rownames(cds@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
colnames(cds@reducedDims$UMAP) <- NULL
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "seurat_clusters",
           label_groups_by_cluster=TRUE,
           label_leaves=TRUE,
           label_branch_points=TRUE) +
    coord_equal()+
    lims(x=c(-11, 10), y=c(-11, 10)) +
    scale_color_manual(values = brewer.pal(7, "Paired"))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) +
    coord_equal()+
    lims(x=c(-11, 10), y=c(-11, 10))

## DEG by pseudotime
deg_pseudotime <- graph_test(cds, neighbor_graph="principal_graph")
pr_deg_ids <- row.names(subset(deg_pseudotime, q_value < 0.05))
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=10^seq(-6,-1))
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$seurat_clusters)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

pheatmap::pheatmap(cds[pr_deg_ids[1:100], order(pData(cds)$seurat_clusters)]@assays[["counts"]], 
                   cluster_cols = F,
                   show_colnames = F,
                   show_rownames = F,
                   scale = "row")

monocle::plot_pseudotime_heatmap(cds[pr_deg_ids,],
                                 cluster_rows = F,
                                 show_rownames = T)

marker_pseudotime <- c("Lef1", "Lmo4", "Rag1", 
                       "Rag2", "Dntt", "Vpreb1", 
                       "Vpreb3", "Ighm", "Ighd",
                       "Cd24a", "Spn", "Kit")
marker_pseudotime_subset <- cds[rowData(cds)$gene_short_name %in% marker_pseudotime, ]
plot_genes_in_pseudotime(marker_pseudotime_subset ,
                         color_cells_by="seurat_clusters",
                         min_expr=0.5,
                         ncol = 3)
