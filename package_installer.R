# this file is responsible for managing all the necessary dependencies and installations required for the analysis. 

# function to validate & install required packages
package_installer <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)  # loads the packages after installing
  }
}


# list of required packages
required_packages <- c("Seurat", "ggplot2", "dplyr", "Rfast2", "SeuratData", "SeuratObject", "patchwork",
                       "plotly", "htmlwidgets", "STdeconvolve", "SpatialExperiment", "cowplot", "Giotto", "FactoMineR", "scran")  # add the packages as needed
package_installer(required_packages)


# installing Giotto
if(!"Giotto" %in% installed.packages()) {
  devtools::install_github("drieslab/Giotto@suite")}


if(!"GiottoData" %in% installed.packages()) {
  devtools::install_github("drieslab/GiottoData")}


genv_exists = checkGiottoEnvironment()  # checking if a Giotto environment already exists
if(!genv_exists){installGiottoEnvironment()}

# x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x
# x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x
# x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x - x