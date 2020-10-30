# Fix PATH to add external dependencies
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/apps/SPAdes/3.13.0/bin/:/apps/HybPiper/:/apps/partitionfinder/2.1.1/", sep = ":"))

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
  verbose = 1,
  parallelism = "clustermq",
  jobs = 20, # Change to match number of cores available!
  cache = skimming_cache,
  seed = 0
)

# - serial
# drake_config(
#   skimming_plan, 
#   verbose = 1, 
#   cache = skimming_cache, 
#   seed = 0
# )
