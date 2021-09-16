## Subset B cells
brain_3prime_bcell <- subset(brain_3prime_singlet, seurat_clusters %in% c(0, 1, 2, 3, 4, 8))
brain_3prime_bcell <- NormalizeData(brain_3prime_bcell)
brain_3prime_bcell <- FindVariableFeatures(brain_3prime_bcell,
                                           selection.method = "vst",
                                           nfeatures = 1500)
brain_3prime_bcell <- ScaleData(brain_3prime_bcell, features = rownames(brain_3prime_bcell))
brain_3prime_bcell <- RunPCA(brain_3prime_bcell,
                             features = VariableFeatures(brain_3prime_bcell))

g1s_list <- c("Cdca7", "Mcm4", "Mcm2", "Rfc2", "Ung", "Mcm6", "Rrm1", "Slbp", 
              "Pcna", "Atad2", "Tipin", "Mcm5", "Uhrf1", "Rpa2", "Dtl", "Prim1", 
              "Fen1", "Hells", "Gmnn", "Pold3", "Nasp", "Chaf1b", "Gins2", "Pola1", 
              "Msh2", "Casp8ap2", "Cdc6", "Ubr7", "Ccne2", "Wdr76", "Tyms", "Cdc45", 
              "Clspn", "Rrm2", "Dscc1", "Rad51", "Usp1", "Exo1", "Blm", "Rad51ap1", 
              "Mlf1ip", "E2f8", "Brip1")
g2m_list <- c("Cbx5", "Aurkb", "Cks1b", "Cks2", "Hn1", "Hmgb2", "Anp32e", "Lbr", 
              "Tmpo", "Top2a", "Tacc3", "Tubb4b", "Ncapd2", "Rangap1", "Cdk1", 
              "Smc4", "Kif20b", "Cdca8", "Ckap2", "Ndc80", "Dlgap5", "Hjurp", 
              "Ckap5", "Bub1", "Ckap2l", "Ect2", "Kif11", "Birc5", "Cdca2", 
              "Nuf2", "Cdca3", "Nusap1", "Ttk", "Aurka", "Mki67", "Fam64a", 
              "Ccnb2", "Tpx2", "Anln", "Kif2c", "Cenpe", "Gtse1", 
              "Kif23", "Cdc20", "Ube2c", "Cenpf", "Cenpa", "Hmmr", "Ctcf", 
              "Psrc1", "Cdc25c", "Nek2", "Gas2l3", "G2e3")
brain_3prime_bcell <- CellCycleScoring(brain_3prime_bcell, 
                                       s.features = g1s_list, 
                                       g2m.features = g2m_list,
                                       set.ident = TRUE)
brain_3prime_bcell <- ScaleData(brain_3prime_bcell, 
                                vars.to.regress = c("S.Score", "G2M.Score"), 
                                features = rownames(brain_3prime_bcell))
brain_3prime_bcell <- RunPCA(brain_3prime_bcell,
                             features = VariableFeatures(brain_3prime_bcell))
ElbowPlot(brain_3prime_bcell, ndims = 30)

brain_3prime_bcell <- FindNeighbors(brain_3prime_bcell, dims = 1:15)
brain_3prime_bcell <- FindClusters(brain_3prime_bcell, resolution = 0.2)
brain_3prime_bcell <- RunUMAP(brain_3prime_bcell, dims = 1:15)

DimPlot(brain_3prime_bcell, label = TRUE) +
    coord_equal()+
    lims(x=c(-11, 9), y=c(-11, 9)) +
    scale_color_manual(values = brewer.pal(7, "Paired"))

DimPlot(brain_3prime_bcell, label = TRUE, split.by = "age") +
    coord_equal()+
    lims(x=c(-11, 9), y=c(-11, 9)) +
    scale_color_manual(values = brewer.pal(7, "Paired"))

## Marker feature plot
p <- FeaturePlot(brain_3prime_bcell, 
                 c("Cd19","Kit", "Spn", "Vpreb1", 
                   "Vpreb3", "Dntt", "Lef1","Rag1", "Rag2",
                   "Top2a", "H2afx","Tuba1b", 
                   "Mki67", "Ccnb2", "Cd24a","Lmo4", 
                   "Ighm", "Ms4a1", 
                   "H2-Aa", "H2-Eb1", "Ighd", 
                   "Bank1", "Prdm1", "Il10"),
                 cols = c("lightgrey", "red"),
                 combine = FALSE,
                 pt.size = 2)
p <- lapply(p, FUN = function(x) AugmentPlot(x) +
                coord_equal()+
                lims(x=c(-11, 9), y=c(-11, 9)))
CombinePlots(plots = p, ncol = 4)

## Marker dot plot
levels(brain_3prime_bcell) <- c(3, 5, 4, 0,1, 2, 6)
DotPlot(brain_3prime_bcell,
        features = rev(c("Cd19","Kit", "Spn", "Vpreb1", 
                         "Vpreb3", "Dntt", "Lef1","Rag1", "Rag2",
                         "Top2a", "H2afx","Tuba1b", 
                         "Mki67", "Ccnb2", "Cd24a","Lmo4", 
                         "Ighm", "Ms4a1", 
                         "H2-Aa", "H2-Eb1", "Ighd", 
                         "Bank1", "Prdm1", "Il10")),
        scale.by = "size") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_color_gradientn(colours = rev(brewer.pal(n = 5, name ="RdYlBu")))

## cluster percentage bar plot
brain_3prime_bcell@meta.data %>%
    mutate(age = ifelse(hash.ID == "hashtag1-TotalSeqB",
                        "2 day", 
                        ifelse(hash.ID %in% c("hashtag2-TotalSeqB", "hashtag3-TotalSeqB"),
                               "2 weeks",
                               "6 weeks"))) %>% 
    group_by(age, seurat_clusters) %>%
    summarise(count = n()) %>%
    ggplot(aes(age, count, fill = seurat_clusters)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = brewer.pal(7, "Paired")) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic()

## DEG
brain_3prime_bcell_marker <- FindAllMarkers(brain_3prime_bcell,
                                              test.use = "MAST",
                                              only.pos = TRUE,
                                              min.pct = 0.25,
                                              logfc.threshold = 0.25)

