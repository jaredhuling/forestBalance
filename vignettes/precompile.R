# Pre-compile vignettes that contain expensive simulations.
# Run this script from the package root before building for CRAN:
#   Rscript vignettes/precompile.R

old_wd <- setwd("vignettes")

vignettes <- c(
  "guide.Rmd.orig",
  "crossfitting.Rmd.orig",
  "augmented.Rmd.orig",
  "performance.Rmd.orig"
)

for (v in vignettes) {
  out <- sub("\\.Rmd\\.orig$", ".Rmd", v)
  message("Knitting ", v, " -> ", out)
  knitr::knit(v, output = out)
}

setwd(old_wd)
