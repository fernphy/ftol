# Setup folders and data for analysis
library(fs)
library(contentid)

# Set up folders ----

dir_create("_targets/user/data_raw/ref_aln")
dir_create("_targets/user/intermediates/blast_sanger")
dir_create("_targets/user/intermediates/iqtree")
dir_create("_targets/user/intermediates/iqtree/plastome")
dir_create("_targets/user/intermediates/iqtree/sanger")
dir_create(paste0("_targets/user/intermediates/iqtree/sanger_", 1:10))
dir_create("_targets/user/intermediates/iqtree/sanger_fast")
dir_create("_targets/user/intermediates/treepl")
dir_create("_targets/user/intermediates/treepl/con")
dir_create("_targets/user/intermediates/treepl/ml")
dir_create("_targets/user/intermediates/treepl/ts")
dir_create("_targets/user/results")
dir_create("ftol_data")

# Fetch data ----

# Download and unzip reference alignments
utils::untar(
  contentid::resolve(
    "hash://sha256/3c37fb9478a8d6d1d5cf12652f01e04c3187db64923be824ca689d924facde18" # nolint
  ),
  exdir = "_targets/user/data_raw/ref_aln"
  )
