source("renv/activate.R")

# Increase time-limit for downloads
options(timeout = max(1000, getOption("timeout")))

options(width = 200)

# Resolve conflicts
conflicted::conflict_prefer("filter", "dplyr", quiet = TRUE)
conflicted::conflict_prefer("select", "dplyr", quiet = TRUE)
conflicted::conflict_prefer("gather", "tidyr", quiet = TRUE)
conflicted::conflict_prefer("map", "purrr", quiet = TRUE)
conflicted::conflict_prefer("resolve", "contentid", quiet = TRUE)

# Return tibble from taxastand functions by default
options(ts_tbl_out = TRUE)