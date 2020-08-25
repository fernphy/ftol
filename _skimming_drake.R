# Load packages, functions, and plan
source("R/packages.R")
source("R/functions.R")
source("R/skimming_plan.R")

# Setup cache
skimming_cache <- new_cache("skimming_cache")
options(rstudio_drake_cache = skimming_cache)

# Specify parallel back-end
options(clustermq.scheduler = "multicore")

# Configure settings for making plan
# (choose either parallel or serial, comment-out the other)

# - parallel
drake_config(
  skimming_plan,
  parallelism = "clustermq",
  jobs = 6, # Change to match number of cores available!
  cache = skimming_cache,
  seed = 0
)

# - serial
# drake_config(plan, verbose = 1, cache = plastid_cache, seed = 0)
