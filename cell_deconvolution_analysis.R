# 2. STDecolvolve
decolvolve_spatial_transcriptomics_analysis <- function(file_path, output_dir, load_data = "filtered", 
                                                        min_lib_size = 100, min_reads = 3, remove_above = 1.0, 
                                                        remove_below = 0.05, n_top_od = 1000, ks_seq = seq(2, 10, by = 1), 
                                                        opt_model = "min", perc_filt = 0.05, beta_scale = 1000) {
  
  outs_dir <- file.path(file_path, "outs")
  filter_data_dir <- file.path(file_path, "outs", "filtered_feature_bc_matrix")
  filter_data <- file.path(file_path, "outs", "filtered_feature_bc_matrix.h5")
  
  # Check if required files and directories exist
  missing_paths <- c()
  for (path in c(outs_dir, filter_data_dir, filter_data)) {
    if (!file.exists(path)) {
      missing_paths <- c(missing_paths, path)
    }
  }
  
  if (length(missing_paths) > 0) {
    stop("Error: The following required files or directories do not exist:\n",
         paste(missing_paths, collapse = "\n"))
  }
  
  # loading the data
  print("Loading the data...")
  se <- SpatialExperiment::read10xVisium(samples = outs_dir, type = "sparse", data = load_data)
  cd <- se@assays@data@listData$counts  # count matrix
  pos <- SpatialExperiment::spatialCoords(se)  # barcodes can be accessed here
  colnames(pos) <- c("y", "x")
  
  counts <- cleanCounts(cd, min.lib.size = min_lib_size, min.reads = min_reads)  # poor genes and barcodes will be removed from cd
  corpus <- restrictCorpus(counts, removeAbove = remove_above, removeBelow = remove_below, nTopOD = n_top_od)  # using top 1000 most significant overdispersed genes by default
  ldas <- fitLDA(t(as.matrix(corpus)), Ks = ks_seq, plot = TRUE, verbose = TRUE)  # fitting LDA model
  optLDA <- optimalModel(models = ldas, opt = opt_model)  # selecting optimal model
  results <- getBetaTheta(optLDA, perc.filt = perc_filt, betaScale = beta_scale)  # Get beta and theta results
  deconProp <- results$theta  # filtering out cell-types
  deconGexp <- results$beta
  
  # Generate plots and save results
  # Plot 1: scatterpie plot of all deconvolved cell types / visualize the barcode proportions of all the deconvolved cell-types in the form of scatterpies
  plt1 <- vizAllTopics(theta = deconProp, pos = pos, r = 35, lwd = 0, showLegend = TRUE, plotTitle = NA) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2)) +
    ggplot2::geom_rect(data = data.frame(pos),
                       ggplot2::aes(xmin = min(x)-90, xmax = max(x)+90,
                                    ymin = min(y)-90, ymax = max(y)+90),
                       fill = NA, color = "black", linetype = "solid", size = 0.5) +
    ggplot2::theme(plot.background = ggplot2::element_blank()) +
    ggplot2::guides(colour = "none")
  
  ggplot2::ggsave(filename = "11_STDeconvolve_deconvolved_cell_types_scatterpie_plot.pdf", plot = plt1, path = output_dir, width = 10, height = 8, dpi = 1000)
  
  # Plot 2: topic-wise scatterpie plots / plotting using vizTopic() for faster plotting by the topic
  ps2 <- lapply(colnames(deconProp), function(celltype) {
    vizTopic(theta = deconProp, pos = pos, topic = celltype, plotTitle = paste0("X", celltype),
             size = 0.9, stroke = 0.2, alpha = 0.5, low = "white", high = "red") +
      ggplot2::guides(colour = "none")
  })
  
  g2 <- cowplot::plot_grid(plotlist = ps2, ncol = 3, nrow = 3, align = "hv")
  ggplot2::ggsave(filename = "12_STDeconvolve_combined_topic_wise_scatterpie_plot.pdf", plot = g2, path = output_dir, width = 12, height = 8, dpi = 1000)
  
  # Plot 3: transcriptional profiles / plotting the top marker genes expressed in the deconvolved cell-types
  geneSymbols <- se@rowRanges@elementMetadata$symbol
  names(geneSymbols) <- names(se@rowRanges)
  colnames(deconGexp) <- geneSymbols[colnames(deconGexp)]
  
  ps3 <- lapply(colnames(deconProp), function(celltype) {
    celltype <- as.numeric(celltype)
    highgexp <- names(which(deconGexp[celltype,] > 3))
    log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
    markers <- names(log2fc)[1]
    
    dat <- data.frame(values = as.vector(log2fc), genes = names(log2fc), order = seq(length(log2fc)))
    dat$selectedLabels <- ""
    dat$selectedLabels[1] <- markers
    
    ggplot2::ggplot(data = dat) +
      ggplot2::geom_col(ggplot2::aes(x = order, y = values, fill = factor(selectedLabels == ""), color = factor(selectedLabels == "")), width = 1) +
      ggplot2::scale_fill_manual(values = c("darkblue", "darkblue")) +
      ggplot2::scale_color_manual(values = c("darkblue", "darkblue")) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(min(log2fc) - 0.3, max(log2fc) + 0.3)) +
      ggplot2::labs(title = paste0("X", celltype), x = "Gene expression rank", y = "log2(FC)") +
      ggplot2::geom_text(ggplot2::aes(x = order+1, y = values-0.1, label = selectedLabels), color = "red") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size=15, color = "black"),
        axis.text.y = ggplot2::element_text(size=15, color = "black"),
        axis.title = ggplot2::element_text(size=15, color = "black"),
        axis.ticks.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size=15),
        panel.grid.major.y = ggplot2::element_line(size = 0.3, colour = "gray80"),
        axis.line = ggplot2::element_line(size = 1, colour = "black"),
        legend.position="none"
      )
  })
  
  g3 <- cowplot::plot_grid(plotlist = ps3, ncol = 3, nrow = 3, align = "hv")
  ggplot2::ggsave(filename = "13_STDeconvolve_combined_transcriptional_profiles_plot.pdf", plot = g3, path = output_dir, width = 12, height = 8, dpi = 1000)
  
  # 3. saving the top genes to .csv
  get_top_genes <- function(deconGexp, deconProp, n = 5) {
    top_genes <- lapply(colnames(deconProp), function(celltype) {
      celltype <- as.numeric(celltype)
      highgexp <- names(which(deconGexp[celltype,] > 3))
      log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
      return(names(log2fc)[1:n])
    })
    return(top_genes)
  }
  
  top_genes_list <- get_top_genes(deconGexp, deconProp)
  top_genes_df <- as.data.frame(do.call(cbind, top_genes_list))
  colnames(top_genes_df) <- paste0("Topic", 1:ncol(top_genes_df))
  rownames(top_genes_df) <- paste0("Gene", 1:nrow(top_genes_df))
  write.csv(top_genes_df, file = file.path(output_dir, "13_STDeconvolve_top_genes_by_celltype.csv"), row.names = TRUE)
  
  # Plot 4: Gene expression visualization / visualizing the expression of some of the genes - vizGeneCounts()
  c <- se@assays@data@listData$counts
  rownames(c) <- geneSymbols[rownames(c)]
  df <- merge(as.data.frame(pos), as.data.frame(t(as.matrix(c))), by = 0)
  
  markerGenes <- unlist(lapply(colnames(deconProp), function(celltype) {
    celltype <- as.numeric(celltype)
    highgexp <- names(which(deconGexp[celltype,] > 3))
    log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
    names(log2fc)[1]
  }))
  
  ps4 <- lapply(markerGenes, function(marker) {
    vizGeneCounts(df = df, gene = marker, size = 0.9, stroke = 0.1, plotTitle = marker, winsorize = 0.05, showLegend = TRUE) +
      ggplot2::guides(colour = "none") +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size=0, color = "black"),
        axis.title = ggplot2::element_text(size=15),
        plot.title = ggplot2::element_text(size=15),
        legend.text = ggplot2::element_text(size = 15, colour = "black"),
        legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 90),
        panel.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(fill = NA, color = "black", size = 2),
        plot.background = ggplot2::element_blank()
      ) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(
        title = "Counts", title.position = "left", title.hjust = 0.5,
        ticks.colour = "black", ticks.linewidth = 2,
        frame.colour = "black", frame.linewidth = 2,
        label.hjust = 0
      ))
  })
  
  g4 <- cowplot::plot_grid(plotlist = ps4, ncol = 3, nrow = 3, align = "hv")
  ggplot2::ggsave(filename = "14_STDeconvolve_combined_genes_viz_geneCounts_plot.pdf", plot = g4, path = output_dir, width = 12, height = 8, dpi = 1000)
  
  return(list(deconProp = deconProp, deconGexp = deconGexp, se = se))  # return
}
