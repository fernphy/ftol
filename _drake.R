# Load packages, functions, and plan
source("R/packages.R")
source("R/functions.R")
source("R/plan.R")

# Setup cache
ftol_cache <- new_cache("ftol_cache")
options(rstudio_drake_cache = ftol_cache)

# Specify parallel back-end
options(clustermq.scheduler = "multicore")

# Configure settings for making plan
# (choose either parallel or serial, comment-out the other)

# - parallel
drake_config(
  plan,
  parallelism = "clustermq",
  jobs = 32, # Change to match number of cores available!
  cache = ftol_cache,
  seed = 0
)

# - serial
# drake_config(plan, verbose = 1, cache = ftol_cache, seed = 0)
