pkgname <- "forestBalance"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('forestBalance')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("compute_balance")
### * compute_balance

flush(stderr()); flush(stdout())

### Name: compute_balance
### Title: Compute covariate balance diagnostics for a set of weights
### Aliases: compute_balance print.forest_balance_diag

### ** Examples





cleanEx()
nameEx("forest_balance")
### * forest_balance

flush(stderr()); flush(stdout())

### Name: forest_balance
### Title: Estimate ATE using forest-based kernel energy balancing
### Aliases: forest_balance

### ** Examples





cleanEx()
nameEx("forest_kernel")
### * forest_kernel

flush(stderr()); flush(stdout())

### Name: forest_kernel
### Title: Compute random forest proximity kernel from a GRF forest
### Aliases: forest_kernel

### ** Examples





cleanEx()
nameEx("get_leaf_node_matrix")
### * get_leaf_node_matrix

flush(stderr()); flush(stdout())

### Name: get_leaf_node_matrix
### Title: Extract leaf node membership matrix from a GRF forest
### Aliases: get_leaf_node_matrix

### ** Examples





cleanEx()
nameEx("kernel_balance")
### * kernel_balance

flush(stderr()); flush(stdout())

### Name: kernel_balance
### Title: Kernel energy balancing weights via closed-form solution
### Aliases: kernel_balance

### ** Examples





cleanEx()
nameEx("leaf_node_kernel")
### * leaf_node_kernel

flush(stderr()); flush(stdout())

### Name: leaf_node_kernel
### Title: Compute random forest proximity kernel from a leaf node matrix
### Aliases: leaf_node_kernel

### ** Examples





cleanEx()
nameEx("simulate_data")
### * simulate_data

flush(stderr()); flush(stdout())

### Name: simulate_data
### Title: Simulate observational study data with confounding
### Aliases: simulate_data

### ** Examples

dat <- simulate_data(n = 500, p = 10, ate = 0)
str(dat)

# True ATE is 0, naive estimate is biased
mean(dat$Y[dat$A == 1]) - mean(dat$Y[dat$A == 0])




cleanEx()
nameEx("summary.forest_balance")
### * summary.forest_balance

flush(stderr()); flush(stdout())

### Name: summary.forest_balance
### Title: Summarize a forest_balance object
### Aliases: summary.forest_balance print.summary.forest_balance

### ** Examples





### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
