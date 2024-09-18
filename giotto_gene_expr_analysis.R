# 3. Giotto
giotto_spatial_transcriptomics_analysis <- function(file_path, output_dir, project_name, xmax_adj=0, xmin_adj=0, ymax_adj=0, ymin_adj=0, giotto_var="Giotto") {
  
  instrs <- createGiottoInstructions(save_dir = output_dir, save_plot = TRUE, show_plot = FALSE)
  
  # checks
  if (!dir.exists(file.path(file_path, "spatial"))) {
    stop("Error: The 'spatial' folder is missing in the specified path.")
  } else {
    if (!file.exists(file.path(file_path, "spatial", "tissue_lowres_image.png"))) {
      stop("Error: The file 'tissue_lowres_image.png' is missing from the spatial folder.")}
    if (!file.exists(file.path(file_path, "spatial", "tissue_positions_list.csv"))) {
      stop("Error: The file 'tissue_positions_list.csv' is missing from the spatial folder.")}
  }
  
  if (!file.exists(file.path(file_path, "raw_feature_bc_matrix.h5"))) {
    stop("Error: The file 'raw_feature_bc_matrix.h5' is missing from the outs folder.")}
  
  # loading the data
  tpp <- read.csv(file.path(file_path, "spatial", "tissue_positions_list.csv"))
  names(tpp) <- NULL
  visium_data <- createGiottoVisiumObject(visium_dir = NULL,
                                          h5_visium_path = file.path(file_path, "raw_feature_bc_matrix.h5"),
                                          h5_tissue_positions_path = file.path(file_path, "spatial", "tissue_positions_list.csv"),
                                          h5_image_png_path = file.path(file_path, "spatial", "tissue_lowres_image.png"),
                                          h5_json_scalefactors_path = file.path(file_path, "spatial", "scalefactors_json.json"),
                                          png_name = "tissue_lowres_image.png",
                                          instructions = instrs)
  
  spatPlot(gobject = visium_data, cell_color = 'in_tissue', show_image = TRUE, point_alpha = 0.7,
           save_param = list(save_name = paste0("15_", giotto_var, "_spatplot_image")))
  
  # updating the giotto plot
  visium_data <- updateGiottoImage(visium_data, image_name = 'image',
                                   xmax_adj = xmax_adj, xmin_adj = xmin_adj,
                                   ymax_adj = ymax_adj, ymin_adj = ymin_adj)
  
  spatPlot(gobject = visium_data, 
           cell_color = 'in_tissue', 
           show_image = TRUE, 
           point_alpha = 0.7,
           save_param = list(save_name = paste0("16_", giotto_var, "_adjusted_spatplot_image")))
  
  spatPlot(gobject = visium_data, cell_color = 'in_tissue', point_size = 2,
           cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
           save_param = list(save_name = paste0("16_", giotto_var, "_adjusted_in_tissue_spatplot_image")))
  
  # subset on spots that were covered by tissue
  metadata <- pDataDT(visium_data)
  in_tissue_barcodes <- metadata[in_tissue == 1]$cell_ID
  visium_data <- subsetGiotto(visium_data, cell_ids = in_tissue_barcodes)
  
  # filtering and nomalization
  visium_data <- filterGiotto(gobject = visium_data, expression_threshold = 1,
                              gene_det_in_min_cells = 10, min_det_genes_per_cell = 500,
                              expression_values = c('raw'), verbose = TRUE)
  visium_data <- normalizeGiotto(gobject = visium_data, scalefactor = 6000, verbose = TRUE)
  visium_data <- addStatistics(gobject = visium_data)
  
  # visuvalize
  spatPlot2D(gobject = visium_data, show_image = TRUE, point_alpha = 0.7,
             save_param = list(save_name = paste0("16_", giotto_var, "_adjusted_spatial_locations")))
  
  spatPlot2D(gobject = visium_data, show_image = TRUE, point_alpha = 0.7,
             cell_color = 'nr_genes', color_as_factor = FALSE,
             save_param = list(save_name = paste0("16_", giotto_var, "_adjusted_num_genes")))
  
  # highly variable genes (HVG) and PCA
  visium_data <- calculateHVG(gobject = visium_data,
                              save_param = list(save_name = paste0("17_", giotto_var, "_HVGplot")))
  gene_metadata <- fDataDT(visium_data)
  featgenes <- gene_metadata[hvg == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$gene_ID
  visium_data <- runPCA(gobject = visium_data, genes_to_use = featgenes,
                        scale_unit = FALSE, center = TRUE, method = "factominer")
  
  screePlot(visium_data, ncp = 30, save_param = list(save_name = paste0("17_", giotto_var, "_screeplot")))
  Giotto::plotPCA(gobject = visium_data, save_param = list(save_name = paste0("18_", giotto_var,"_PCA_reduction")))
  
  # UMAP and tSNE
  visium_data <- Giotto::runUMAP(visium_data, dimensions_to_use = 1:10)
  Giotto::plotUMAP(gobject = visium_data, save_param = list(save_name = paste0("18_", giotto_var, "_UMAP_reduction")))
  
  visium_data <- Giotto::runtSNE(visium_data, dimensions_to_use = 1:10)
  Giotto::plotTSNE(gobject = visium_data, save_param = list(save_name = paste0("18_", giotto_var, "_tSNE_reduction")))
  
  # Clustering
  visium_data <- createNearestNetwork(gobject = visium_data, dimensions_to_use = 1:10, k = 15)  # snn network by default
  visium_data <- doLeidenCluster(gobject = visium_data, resolution = 0.4, n_iterations = 1000)  # Leiden clustering
  
  plotUMAP(gobject = visium_data, cell_color = 'leiden_clus', show_NN_network = TRUE, point_size = 2.5,
           save_param = list(save_name = paste0("19_", giotto_var, "_UMAP_leiden")))
  
  # spatial expression and dimensional plotting
  spatDimPlot(gobject = visium_data, cell_color = 'leiden_clus',
              dim_point_size = 2, spat_point_size = 2.5,
              save_param = list(save_name = paste0("19_", giotto_var, "_covis_leiden")))
  
  spatDimPlot(gobject = visium_data, cell_color = 'nr_genes', color_as_factor = FALSE,
              dim_point_size = 2, spat_point_size = 2.5,
              save_param = list(save_name = paste0("19_", giotto_var, "_num_genes")))
  
  # visualize UMAP and spatial results
  spatDimPlot(gobject = visium_data, cell_color = 'leiden_clus', spat_point_shape = 'voronoi',
              save_param = list(save_name = paste0("19_", giotto_var, "_covis_leiden_voronoi")))
  
  # heatmap and dendrogram
  showClusterHeatmap(gobject = visium_data, cluster_column = 'leiden_clus',
                     save_param = list(save_name = paste0("20_", giotto_var, "_num_genes")))
  showClusterDendrogram(visium_data, h = 0.5, rotate = T, cluster_column = 'leiden_clus',
                        save_param = list(save_name = paste0("20_", giotto_var, "_num_genes")))
  
  #gini
  gini_markers_subclusters = findMarkers_one_vs_all(gobject = visium_data,
                                                    method = 'gini',
                                                    expression_values = 'normalized',
                                                    cluster_column = 'leiden_clus',
                                                    min_genes = 20,
                                                    min_expr_gini_score = 0.5,
                                                    min_det_gini_score = 0.5)
  topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$genes
  write.csv(gini_markers_subclusters, file.path(output_dir, "21_gini_markers_subclusters.csv"))
  
  #####################
  # violinplot
  violinPlot(visium_data, genes = unique(topgenes_gini), cluster_column = 'leiden_clus',
             strip_text = 8, strip_position = 'right',
             save_param = list(save_name = paste0("21_", giotto_var, '_violinplot_gini'), base_width = 5, base_height = 10))
  
  
  # cluster heatmap
  plotMetaDataHeatmap(visium_data, selected_genes = topgenes_gini,
                      metadata_cols = c('leiden_clus'), x_text_size = 10, y_text_size = 10,
                      save_param = list(save_name = paste0("21_", giotto_var, "_metaheatmap_gini")))
  
  # umap plots
  dimGenePlot2D(visium_data, expression_values = 'scaled',
                genes = gini_markers_subclusters[, head(.SD, 1), by = 'cluster']$genes,
                cow_n_col = 3, point_size = 1,
                save_param = list(save_name = paste0("21_", giotto_var, "_gini_umap"), base_width = 8, base_height = 5))
  
  scran_markers_subclusters = findMarkers_one_vs_all(gobject = visium_data, method = 'scran',
                                                     expression_values = 'normalized', cluster_column = 'leiden_clus')
  topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$genes
  
  write.csv(scran_markers_subclusters, file.path(output_dir, "22_scran_markers_subclusters.csv"))
  
  
  # violinplot
  violinPlot(visium_data, genes = unique(topgenes_scran), cluster_column = 'leiden_clus',
             strip_text = 10, strip_position = 'right',
             save_param = list(save_name = paste0("22_", giotto_var, "_violinplot_scran"), base_width = 5))
  
  # cluster heatmap
  plotMetaDataHeatmap(visium_data, selected_genes = topgenes_scran, metadata_cols = c('leiden_clus'),
                      save_param = list(save_name = paste0("22_", giotto_var, "_metadata_heatmap")))
  
  # umap plots
  dimGenePlot(visium_data, expression_values = 'scaled',
              genes = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$genes,
              cow_n_col = 3, point_size = 1,
              save_param = list(save_name = paste0("22_", giotto_var, "_scran_umap")), base_width = 8, base_height = 5)
  
  # Save the processed object as an RDS file
  saveRDS(visium_data, file = file.path(output_dir, paste0("22_", giotto_var, "_", project_name, ".rds")))
  
  return(visium_data)
}

# x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x
# x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x
# x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x
