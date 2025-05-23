# loading custom functions
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### NOTE:
# /path/to/sample/outs is where my 10x-visium data is present:
# Seurat: /path/to/sample/outs
# STDconvolve: /path/to/sample
# Giotto: /path/to/sample/outs

files <- list.files()
cat("Files in the config dir:\n", paste(files, collapse = "\n"))

print(" ")
print(" ")
print("Loading necessary files:")


source("package_installer.R")
source("seurat_analysis.R")
source("cell_deconvolution_analysis.R")
source("giotto_gene_expr_analysis.R")


# Usage:
automated_seurat_analysis_func(
  data_dir = "/path/to/sample/outs",
  project_name = "CUSTOM_PROJECT_NAME",
  output_dir = "/path/to/the/output_directory"
)

automated_seurat_spatial_analysis_func(
  data_dir = "/path/to/sample/outs",
  file_name = "filtered_feature_bc_matrix.h5", 
  slice_no = "SLIDE_NO", 
  output_dir = "/path/to/the/output_directory"
)

decolvolve_spatial_transcriptomics_analysis(file_path = "/path/to/sample",
                                            output_dir = "/path/to/the/output_directory")

giotto_spatial_transcriptomics_analysis(file_path = "/path/to/sample/outs",
                                        output_dir = "/path/to/sample",
                                        project_name = "CUSTOM_PROJECT_NAME",
                                        xmax_adj = 1300, xmin_adj = 1400, ymax_adj = 1100, ymin_adj = 1200)
# xmax_adj, xmin_adj, ymax_adj, ymin_adj vary based on the tissue, however we used the following values respectively;
# xmax_adj = 1300, xmin_adj = 1400, ymax_adj = 1100, ymin_adj = 1200

# x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x
# x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x
# x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x
