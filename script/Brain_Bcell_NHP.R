brain_bcell <- subset(brain_mmul_merge,
                      subset = seurat_clusters %in% c(0, 1, 2, 3, 5, 6, 11))
brain_bcell <- NormalizeData(brain_bcell)
brain_bcell <- FindVariableFeatures(brain_bcell,
                                    selection.method = "vst",
                                    nfeatures = 2000)
brain_bcell <- ScaleData(brain_bcell, features = rownames(brain_bcell))
brain_bcell <- RunPCA(brain_bcell,
                      features = VariableFeatures(brain_bcell))

g1s_list <- read.csv("signature_gene/g1s_mmul.csv", stringsAsFactors = FALSE)
g1s_list <- g1s_list$Gene
g2m_list <- read.csv("signature_gene/g2m_mmul.csv", stringsAsFactors = FALSE)
g2m_list <- g2m_list$Gene

brain_bcell <- CellCycleScoring(brain_bcell, 
                                s.features = g1s_list, 
                                g2m.features = g2m_list,
                                set.ident = TRUE)
brain_bcell <- ScaleData(brain_bcell, 
                         vars.to.regress = c("S.Score", "G2M.Score"), 
                         features = rownames(brain_bcell))
brain_bcell <- RunPCA(brain_bcell,
                      features = VariableFeatures(brain_bcell))
ElbowPlot(brain_bcell, ndims = 50)

brain_bcell <- FindNeighbors(brain_bcell, dims = 1:30)
brain_bcell <- FindClusters(brain_bcell, resolution = 0.5)
brain_bcell <- RunUMAP(brain_bcell, dims = 1:30)
DimPlot(brain_bcell, label = TRUE)


## custom dot plot
PercentAbove <- function(x, threshold) {
    return(length(x = x[x > threshold]) / length(x = x))
}
dotplot_custom <- function (object, assay = NULL, features, cols = c("lightgrey", 
                                                                     "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, 
                            dot.scale = 6, group.by = NULL, split.by = NULL, scale.by = "radius", 
                            scale.min = NA, scale.max = NA) 
{
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    scale.func <- switch(EXPR = scale.by, size = scale_size, 
                         radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    data.features <- FetchData(object = object, vars = features)
    data.features$RAG1 <- 0
    data.features$RAG2 <- 0
    data.features$DNTT <- 0
    data.features$id <- if (is.null(x = group.by)) {
        Idents(object = object)
    }
    else {
        object[[group.by, drop = TRUE]]
    }
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]]
        if (length(x = unique(x = splits)) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                            "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident, 
                                  1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
            return(mean(x = expm1(x = x)))
        })
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                         threshold = 0)
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                             FUN = function(x) {
                                 data.use <- data.plot[data.plot$features.plot == 
                                                           x, "avg.exp"]
                                 data.use <- scale(x = data.use)
                                 data.use <- MinMax(data = data.use, min = col.min, 
                                                    max = col.max)
                                 return(data.use)
                             })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                             breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                      levels = rev(x = features))
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
        splits.use <- vapply(X = strsplit(x = as.character(x = data.plot$id), 
                                          split = "_"), FUN = "[[", FUN.VALUE = character(length = 1L), 
                             2)
        data.plot$colors <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("grey", 
                                               color))(20)[value])
        }, color = cols[splits.use], value = avg.exp.scaled)
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled", 
                       no = "colors")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                          y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
                                                                                                       color = color.by)) + scale.func(range = c(0, dot.scale), 
                                                                                                                                       limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                                 axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
        labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
                                        yes = "Identity", no = "Split Identity")) + 
        theme_cowplot()
    if (!is.null(x = split.by)) {
        plot <- plot + scale_color_identity()
    }
    else if (length(x = cols) == 1) {
        plot <- plot + scale_color_distiller(palette = cols)
    }
    else {
        plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    }
    if (is.null(x = split.by)) {
        plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
    }
    return(plot)
}

## Marker dot plot
dotplot_custom(brain_bcell,
               features = rev(toupper(c("Cd19","Kit", "Spn", "Vpreb1", 
                                        "Vpreb3", "Dntt", "Lef1", "Rag1", "Rag2",
                                        "Top2a", "H2afx","Tuba1b", 
                                        "Mki67", "Ccnb2", "Lmo4", 
                                        "Ighm", "CD83",  
                                        "Ms4a1","S100A4", "GPR183", "CD86", "ITGAX", "CD44", "FCER2", 
                                        "CD52", "CD1C", "Bank1", "CR2", "IL16", "BCL6", "CD22","Prdm1", "JCHAIN", 
                                        "MZB1", "XBP1", "Il10" ))),
               scale.by = "size") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_color_gradientn(colours = rev(brewer.pal(n = 5, name ="RdYlBu")))


## Marker
brain_bcell_markers <- FindAllMarkers(brain_bcell,
                                      test.use = "MAST", 
                                      latent.vars = "orig.ident",
                                      only.pos = TRUE, 
                                      min.pct = 0.25, 
                                      logfc.threshold = 0.25)

## DEG between each pair of clusters
deg_pairwise <- apply(combn(levels(brain_bcell$seurat_clusters), 2), 
                      2, 
                      function(x) FindMarkers(brain_bcell,
                                              ident.1 = x[1],
                                              ident.2 = x[2],
                                              test.use = "MAST",
                                              latent.vars = "orig.ident",
                                              min.pct = 0.25,
                                              logfc.threshold = 0.25))
names(deg_pairwise) <- apply(combn(levels(brain_bcell$seurat_clusters), 2), 
                             2, 
                             function(x) paste0("cluster_", x[1],"_vs_",x[2], "_wo_cc"))
lapply(names(deg_pairwise), function(df) write.csv(deg_pairwise[[df]], file=paste0("deg_", df, ".csv")))

FeaturePlot(brain_bcell, 
            c("LEF1", "LMO4", "RAG1", 
              "RAG2", "DNTT", "VPREB1", 
              "VPREB3", "IGHM", "MKI67",
              "SPN", "KIT", "PRDM1"), 
            ncol = 3)
