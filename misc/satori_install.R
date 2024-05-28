.libPaths("/nobackup/users/ydeng9/packages/R")
install.packages("SeuratObject", repos = "https://cran.r-project.org")
install.packages("uwot", repos = "https://cran.r-project.org")

#install.packages("Matrix", repos="http://R-Forge.R-project.org")
#install.packages("https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-51.6.tar.gz", repos = NULL, type = "source")
#install.packages("https://cran.r-project.org/src/contrib/Archive/mgcv/mgcv_1.8-31.tar.gz", repos = NULL, type = "source")
#install.packages("https://cran.r-project.org/src/contrib/Archive/digest/digest_0.6.18.tar.gz", repos = NULL, type = "source")
#install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.3.3.tar.gz", repos = NULL, type = "source")
#install.packages("remotes", repos = "https://cran.r-project.org")
remotes::install_version("uwot", version = "0.1.10")
remotes::install_version("SeuratObject", version = "4.0.0")
remotes::install_version("Seurat", version = "3.2.3")

#install.packages("Seurat", repos = "https://cran.r-project.org")

