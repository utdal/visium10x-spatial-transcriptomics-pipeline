# 1. Seurat Analysis
automated_seurat_analysis_func <- function(data_dir, project_name, output_dir) {
  setwd(data_dir)  # setting working directory
  drg_spatial.data <- Read10X(data.dir = file.path(data_dir, "filtered_feature_bc_matrix"))  # loading the data
  drg_spatial <- CreateSeuratObject(counts = drg_spatial.data, project = "drgspatial", min.cells = 3, min.features = 200)  # creating a seurat object
  
  if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}
  
  pdf(file = file.path(output_dir, "1_seurat_ViolinPlot.pdf"), width = 12, height = 10)
  print(VlnPlot(drg_spatial, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()
  
  drg_spatial[["percent.mt"]] <- PercentageFeatureSet(drg_spatial, pattern = "^MT-")  # calculating MT %
  pdf(file = file.path(output_dir, "1_seurat_ViolinPlot_MT_percent.pdf"), width = 12, height = 10)
  print(VlnPlot(drg_spatial, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()
  
  pdf(file = file.path(output_dir, "1_seurat_Feature_scatter_plot.pdf"), width = 12, height = 10)
  plot1 <- FeatureScatter(drg_spatial, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(drg_spatial, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
  dev.off()
  
  drg_spatial <- NormalizeData(drg_spatial)  # normalization
  
  # identifying variable features
  drg_spatial <- FindVariableFeatures(drg_spatial, selection.method = "vst")
  
  # plotting variable features
  pdf(file = file.path(output_dir, "2_seurat_FindVariableFeatures_plot.pdf"), width = 12, height = 10)
  top10 <- head(VariableFeatures(drg_spatial), 10)
  plot1 <- VariableFeaturePlot(drg_spatial)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  print(plot1 + plot2)
  dev.off()
  
  # scaling the data
  all.genes <- rownames(drg_spatial)
  drg_spatial <- ScaleData(drg_spatial, features = all.genes)
  
  # running PCA
  drg_spatial <- RunPCA(drg_spatial, features = VariableFeatures(object = drg_spatial))
  pdf(file = file.path(output_dir, "3_seurat_PCA_VizDimLoadings_plot.pdf"), width = 12, height = 10)
  print(VizDimLoadings(drg_spatial, dims = 1:4, reduction = "pca"))
  dev.off()
  
  pdf(file = file.path(output_dir, "4_seurat_Dim_plot.pdf"), width = 12, height = 10)
  print(DimPlot(drg_spatial, reduction = "pca") + NoLegend())
  dev.off()
  
  pdf(file = file.path(output_dir, "4_seurat_Dim_Heatmap_plot.pdf"), width = 12, height = 10)
  print(DimHeatmap(drg_spatial, dims = 1:6, cells = 500, balanced = TRUE))
  dev.off()
  
  pdf(file = file.path(output_dir, "4_seurat_Elbow_plot.pdf"), width = 12, height = 10)
  print(ElbowPlot(drg_spatial))
  dev.off()
  
  # clustering cells
  drg_spatial <- FindNeighbors(drg_spatial, dims = 1:10)
  drg_spatial <- FindClusters(drg_spatial, resolution = 0.5)
  
  # running non-linear dim reduction (UMAP/tSNE)
  drg_spatial <- RunUMAP(drg_spatial, dims = 1:10)
  DimPlot(drg_spatial, reduction = "umap")
  pdf(file = file.path(output_dir, "5_seurat_Dim_Clustering_UMAP_plot.pdf"), width = 12, height = 10)
  dim_plot <- DimPlot(drg_spatial, reduction = "umap")
  print(dim_plot)
  dev.off()
  
  # saving the file post clustering/UMAP
  saveRDS(drg_spatial, file = file.path(output_dir, paste0("5_seurat_", project_name, ".rds")))
  
  # finding DE expressed features / cluster biomarkers
  drg_spatial.markers <- FindAllMarkers(drg_spatial, only.pos = TRUE)
  drg_spatial.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)
  
  top2_markers <- drg_spatial.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)  # Get the top 2 genes per cluster based on avg_log2FC
  
  write.csv(top2_markers, file.path(output_dir, paste0("5_seurat_", project_name, "_marker_genes_ordered_by_log2FC.csv")))
  
  pdf(file = file.path(output_dir, "5_seurat_top2_FindAllMarkerGenes_from_each_Cluster.pdf"), width = 12, height = 10)
  vln_plot  <- VlnPlot(drg_spatial, features=top2_markers$gene)
  print(vln_plot)
  dev.off()
}


automated_seurat_spatial_analysis_func <- function(data_dir, file_name, slice_no, output_dir) {
  setwd(data_dir)  # Setting working directory
  
  # Load spatial data and create Seurat object
  data <- Load10X_Spatial(data.dir = data_dir, filename = file_name, 
                          assay = "Spatial", slice = slice_no,
                          bin.size = NULL, filter.matrix = TRUE, 
                          to.upper = FALSE, image = NULL)
  
  # Violin plot                      
  pdf(file = file.path(output_dir, "6_spatial_ViolinPlot.pdf"), width = 12, height = 10)
  print(VlnPlot(data, features = "nCount_Spatial", pt.size = 0.1) + NoLegend())
  dev.off()
  
  # Spatial feature plot             
  pdf(file = file.path(output_dir, "6_spatial_SpatialFeaturePlot.pdf"), width = 12, height = 10)
  print(SpatialFeaturePlot(data, features = "nCount_Spatial") + theme(legend.position = "right"))
  dev.off()
  
  pdf(file = file.path(output_dir, "6_spatial_nCount_Vln_and_SpatialFeaPlot.pdf"), width = 17, height = 12)
  p1 <- VlnPlot(data, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  p2 <- SpatialFeaturePlot(data, features = "nCount_Spatial") + theme(legend.position = "right")
  print(wrap_plots(p1, p2))
  dev.off()
  
  # Apply SCTransform
  data <- SCTransform(data, assay = "Spatial", verbose = FALSE)
  
  # Dim reduction and clustering
  data <- RunPCA(data, assay = "SCT", verbose = FALSE)
  data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
  data <- FindClusters(data, verbose = FALSE)
  data <- RunUMAP(data, reduction = "pca", dims = 1:30)
  
  # UMAP and spatial UMAP plot
  pdf(file = file.path(output_dir, "7_spatial_UMAP_clustering_plot.pdf"), width = 18, height = 10)
  p1 <- DimPlot(data, reduction = "umap", label = TRUE)
  p2 <- SpatialDimPlot(data, label = TRUE, label.size = 3)
  print(p1 + p2)
  dev.off()
  
  # Spatial dim plot by cluster
  pdf(file = file.path(output_dir, "8_spatial_SpatialDimPlot_by_cluster.pdf"), width = 12, height = 10)
  clusters <- Idents(data)
  num_clusters <- length(unique(clusters))
  cluster_indices <- as.character(seq(0, num_clusters - 1))
  print(SpatialDimPlot(data, cells.highlight = CellsByIdentities(object = data, idents = cluster_indices), facet.highlight = TRUE, ncol = 3))
  dev.off()
  
  # interactive plots
  # interactive_spatial_dim_plot <- SpatialDimPlot(data, interactive = TRUE)
  # interactive_feature_plot <- SpatialFeaturePlot(data, features = "Ttr", interactive = TRUE)
  # linked_dim_plot <- LinkedDimPlot(data)
  
  # Save interactive plots to HTML files
  # saveWidget(interactive_spatial_dim_plot, file = file.path(output_dir, "interactive_spatial_dim_plot.html"))
  # saveWidget(interactive_feature_plot, file = file.path(output_dir, "interactive_feature_plot.html"))
  # saveWidget(linked_dim_plot, file = file.path(output_dir, "linked_dim_plot.html"))
  
  # DE analysis
  de_markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(de_markers, file.path(output_dir, "9_spatial_DE_Markers.csv"))
  top2_genes_per_cluster <- de_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
  write.csv(top2_genes_per_cluster, file.path(output_dir, "9_top2_genes_per_cluster_based_on_Log2FC_spatial_FindAllMarkers.csv"))
  pdf(file = file.path(output_dir, "10_spatial_DE_Markers_SpatialFeaturePlot.pdf"), width = 15, height = 30)
  print(SpatialFeaturePlot(data,features = top2_genes_per_cluster$gene, alpha = c(0.1, 1), ncol = 2))
  dev.off()
}

# x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x
# x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x
# x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x