# Pre-compile the performance vignette (which takes ~5 minutes to run)
# Run this script from the package root before building for CRAN:
#   Rscript vignettes/precompile.R

old_wd <- setwd("vignettes")
knitr::knit("performance.Rmd.orig", output = "performance.Rmd")
setwd(old_wd)
