library(Seurat)
library(Matrix)
library(tidyverse)
library(Matrix.utils)
library(pheatmap)
library(multcomp)
library(CountClust)
library(ggbeeswarm)
library(RColorBrewer)

read_10x_matrix <- function(matrix_dir, project_id){
    barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
    features.path <- paste0(matrix_dir, "features.tsv.gz")
    matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
    mat <- readMM(file = matrix.path)
    feature.names = read.delim(features.path, 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path, 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    colnames(mat) = barcode.names$V1
    rownames(mat) = feature.names$V2
    ## some ensembl gene id correspond to same gene symbol
    ## thus we need to combine these rows 
    mat <- aggregate.Matrix(mat, rownames(mat), FUN="sum")
    seurat_object <- CreateSeuratObject(counts = mat, 
                                        project = project_id, 
                                        min.cells = 3, 
                                        min.features = 200)
    return(seurat_object)
}

## Read count matrix
brain_mmul_1 <- read_10x_matrix("data/filtered_feature_bc_matrix_DuB1/", "brain_mmul_1")
brain_mmul_2 <- read_10x_matrix("data/filtered_feature_bc_matrix_DuB2/", "brain_mmul_2")
brain_mmul_3 <- read_10x_matrix("data/filtered_feature_bc_matrix_DuB3/", "brain_mmul_3")

## Merge seurat object
cell_prefix_brain_mmul <- c("brain_mmul_1", "brain_mmul_2", "brain_mmul_3")
brain_mmul_merge <- merge(brain_mmul_1, 
                          c(brain_mmul_2, brain_mmul_3),
                          add.cell.ids = cell_prefix_brain_mmul)
MT_gene <- c("ATP6","ATP8","CYTB","COX1","COX2","COX3","ND1","ND2","ND3","ND4","ND4L","ND5","ND6")
MT_filt <- MT_gene %in% rownames(brain_mmul_merge[['RNA']]@counts)
brain_mmul_merge[["percent.mt"]] <- PercentageFeatureSet(brain_mmul_merge, 
                                                         features = MT_gene[MT_filt])

## Check cutoff for count and feature
## nCount_RNA
brain_mmul_merge@meta.data %>% 
    filter(nCount_RNA<30000) %>% 
    ggplot(aes(nCount_RNA)) +
    geom_histogram(binwidth = 100, boundary = 0) +
    facet_wrap(~orig.ident) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 30)) +
    geom_vline(xintercept = c(1500, 10000)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))

## nFeature_RNA
brain_mmul_merge@meta.data %>% 
    filter(nFeature_RNA<5000) %>% 
    ggplot(aes(nFeature_RNA)) +
    geom_histogram(binwidth = 20, boundary = 0) +
    facet_wrap(~orig.ident) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 30)) +
    geom_vline(xintercept = c(600)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))

## percent.mt
brain_mmul_merge@meta.data %>% 
    ggplot(aes(percent.mt)) +
    geom_histogram(binwidth = 0.01, boundary = 0) +
    facet_wrap(~orig.ident) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 30)) +
    geom_vline(xintercept = 5) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))

## Subset
brain_mmul_merge <- subset(brain_mmul_merge,
                           subset = nFeature_RNA > 600 & 
                               nCount_RNA < 12000 & 
                               nCount_RNA > 1500 &
                               percent.mt < 5)

brain_mmul_merge <- NormalizeData(brain_mmul_merge)
brain_mmul_merge <- FindVariableFeatures(brain_mmul_merge, 
                                         selection.method = "vst", 
                                         nfeatures = 2000)
all.genes_brain_mmul_merge <- rownames(brain_mmul_merge)
brain_mmul_merge  <- ScaleData(brain_mmul_merge, features = all.genes_brain_mmul_merge)

brain_mmul_merge <- RunPCA(brain_mmul_merge, features = VariableFeatures(object = brain_mmul_merge))
ElbowPlot(brain_mmul_merge, ndims = 50)

brain_mmul_merge <- FindNeighbors(brain_mmul_merge, dims = 1:25)
brain_mmul_merge <- FindClusters(brain_mmul_merge, resolution = 0.4)
brain_mmul_merge <- RunUMAP(brain_mmul_merge, dims = 1:25)
DimPlot(brain_mmul_merge, label = TRUE) 

brain_mmul_merge_markers <- FindAllMarkers(brain_mmul_merge,
                                           test.use = "MAST", 
                                           latent.vars = "orig.ident",
                                           only.pos = TRUE)



FeaturePlot(brain_mmul_merge, 
            c("LEF1", "LMO4", "RAG1", 
              "RAG2", "DNTT", "VPREB1", 
              "VPREB3", "IGHM", "MKI67",
              "SPN", "KIT", "PRDM1"), 
            ncol = 3)

DotPlot(subset(brain_mmul_merge, subset = seurat_clusters %in% c(0, 1, 2, 3, 5, 6, 8, 11)),
        features = rev(toupper(c("Cd19","Kit", "Spn", "Vpreb1", 
                                 "Vpreb3", "Dntt", "Lef1", "Rag1", "Rag2",
                                 "Top2a", "H2afx","Tuba1b", 
                                 "Mki67", "Ccnb2", "Lmo4", 
                                 "Ighm", "CD83",  
                                 "Ms4a1","S100A4", "GPR183", "CD86", "ITGAX", "CD44", "FCER2", 
                                 "CD52", "CD1C", "Bank1", "CR2", "IL16","BCL6", "CD22","Prdm1", "JCHAIN", 
                                 "MZB1", "XBP1", "Il10"))),
        scale.by = "size") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_color_gradientn(colours = rev(brewer.pal(n = 5, name ="RdYlBu")))
