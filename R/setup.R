# Setup folders and data for analysis
library(fs)

# Set up folders ----

dir_create("_targets/user/data_raw")
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
dir_create("_targets/user/results/ftolr")

# Download data from FigShare ----

# Download input data to a temporary zip file
temp_file <- tempfile(fileext = ".zip")
download.file(
  "https://figshare.com/ndownloader/articles/19474316/versions/1",
  destfile = temp_file)

# Unzip to targets/user/data_raw
utils::unzip(temp_file, exdir = "_targets/user/data_raw")

# Untar a tar achive that was inside the zipped file
utils::untar(
  "_targets/user/data_raw/ref_aln.tar.gz",
  exdir = "_targets/user/data_raw/ref_aln")

file_delete(temp_file)
