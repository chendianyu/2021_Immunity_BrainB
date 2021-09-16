library(Seurat)
library(Matrix)
library(tidyverse)
library(Matrix.utils)
library(pheatmap)
library(multcomp)
library(CountClust)
library(RColorBrewer)
library(ggbeeswarm)

read_10x_matrix <- function(matrix_dir){
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
    return(mat)
}

brain_3prime_mat <- read_10x_matrix("data/filtered_feature_bc_matrix/")
brain_3prime_hashtag <- brain_3prime_mat[str_detect(rownames(brain_3prime_mat), "hashtag"), ]
brain_3prime_umis <- brain_3prime_mat[!str_detect(rownames(brain_3prime_mat), "hashtag"), ]

# HTO demultiplexing
brain_3prime <- CreateSeuratObject(counts = brain_3prime_umis,
                                   project = "brain_3prime")
brain_3prime[['percent.mt']] <- PercentageFeatureSet(brain_3prime, pattern = "^mt-")
brain_3prime <- NormalizeData(brain_3prime)
brain_3prime <- FindVariableFeatures(brain_3prime, 
                                     selection.method = "mean.var.plot")
brain_3prime <- ScaleData(brain_3prime, features = VariableFeatures(brain_3prime))

## Adding HTO data as an independent assay
brain_3prime[['HTO']] <- CreateAssayObject(counts = brain_3prime_hashtag)
brain_3prime <- NormalizeData(brain_3prime, assay = "HTO", normalization.method = "CLR")
brain_3prime <- HTODemux(brain_3prime, assay = "HTO", positive.quantile = 0.99)

## Generate a tSNE embedding for HTOs
## Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = brain_3prime, assay = "HTO"))))
# Calculate tSNE embeddings with a distance matrix
brain_3prime <- RunTSNE(brain_3prime, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(brain_3prime) 

# Filter out doublets and negative cells
brain_3prime_singlet <- subset(brain_3prime, HTO_classification.global == "Singlet")
brain_3prime_singlet <- brain_3prime_singlet[rowSums(as.matrix(brain_3prime_singlet[["RNA"]]@counts)>0)>3, ]

## QC
## nCount_RNA
brain_3prime_singlet@meta.data %>% 
    filter(nCount_RNA<60000) %>% 
    ggplot(aes(nCount_RNA)) +
    geom_histogram(binwidth = 200, boundary = 0) +
    facet_wrap(~hash.ID) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 30)) +
    geom_vline(xintercept = c(20000)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))

## nFeature_RNA
brain_3prime_singlet@meta.data %>% 
    ggplot(aes(nFeature_RNA)) +
    geom_histogram(binwidth = 20, boundary = 0) +
    facet_wrap(~hash.ID) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 30)) +
    geom_vline(xintercept = c(1200)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))

## percent.mt
brain_3prime_singlet@meta.data %>% 
    filter(percent.mt < 20) %>% 
    ggplot(aes(percent.mt)) +
    geom_histogram(binwidth = 0.1, boundary = 0) +
    facet_wrap(~hash.ID) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 30)) +
    geom_vline(xintercept = c(10)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))

## Subset
brain_3prime_singlet <- subset(brain_3prime_singlet,
                           subset = nFeature_RNA > 1200 & nCount_RNA < 20000 & percent.mt < 10)
## Add age information
brain_3prime_singlet$age <- ifelse(brain_3prime_singlet$hash.ID == "hashtag1-TotalSeqB",
                                   "2 day", 
                                   ifelse(brain_3prime_singlet$hash.ID %in% c("hashtag2-TotalSeqB", "hashtag3-TotalSeqB"),
                                          "2 weeks",
                                          "6 weeks"))

brain_3prime_singlet <- NormalizeData(brain_3prime_singlet)
brain_3prime_singlet <- FindVariableFeatures(brain_3prime_singlet, 
                                         selection.method = "vst", 
                                         nfeatures = 2000)
all.genes_brain_3prime_singlet <- rownames(brain_3prime_singlet)
brain_3prime_singlet  <- ScaleData(brain_3prime_singlet, features = all.genes_brain_3prime_singlet)

brain_3prime_singlet <- RunPCA(brain_3prime_singlet, features = VariableFeatures(object = brain_3prime_singlet))
ElbowPlot(brain_3prime_singlet, ndims = 50)

brain_3prime_singlet <- FindNeighbors(brain_3prime_singlet, dims = 1:20)
brain_3prime_singlet <- FindClusters(brain_3prime_singlet, resolution = 0.3)
brain_3prime_singlet <- RunUMAP(brain_3prime_singlet, dims = 1:20)
DimPlot(brain_3prime_singlet, label = TRUE) +
    coord_equal()+
    lims(x=c(-13, 16), y=c(-13, 16)) +
    scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(15))

AugmentPlot(DimPlot(brain_3prime_singlet, group.by = "orig.ident", cols = "lightgrey")) +
    coord_equal()+
    lims(x=c(-13, 16), y=c(-13, 16))
AugmentPlot(FeaturePlot(brain_3prime_singlet, "Cd19", cols = c("lightgrey", "red"), pt.size = 2)) +
    coord_equal()+
    lims(x=c(-13, 16), y=c(-13, 16))

p <- DimPlot(brain_3prime_singlet, label = TRUE, split.by = "hash.ID", pt.size = 1, combine = FALSE) 
p <- lapply(X = p, FUN = AugmentPlot)
CombinePlots(plots = p, ncol = 5)

brain_3prime_singlet_marker <- FindAllMarkers(brain_3prime_singlet,
                                              test.use = "MAST",
                                              only.pos = TRUE,
                                              min.pct = 0.25,
                                              logfc.threshold = 0.25)
