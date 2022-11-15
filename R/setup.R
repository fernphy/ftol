# Setup folders and data for analysis
library(fs)

# Set up folders ----

dir_create("_targets/user/data_raw/ref_aln")
dir_create("_targets/user/data_raw/restez/sql_db")
dir_create("_targets/user/intermediates/blast_sanger")
dir_create("_targets/user/intermediates/iqtree")
dir_create("_targets/user/intermediates/ref_seqs")
dir_create("_targets/user/intermediates/iqtree/plastome")
dir_create(paste0("_targets/user/intermediates/iqtree/sanger_", 1:10))
dir_create("_targets/user/intermediates/iqtree/sanger_fast")
dir_create("_targets/user/intermediates/treepl")
dir_create("_targets/user/intermediates/treepl/con")
dir_create("_targets/user/intermediates/treepl/ml")
dir_create("_targets/user/intermediates/treepl/ts")
dir_create("_targets/user/results")
dir_create("ftol_data")

# Fetch data ----

# Download and unzip reference alignments from figshare
archive::archive_extract(
  contentid::resolve(
    "hash://sha256/388b53201a8626d4b41851e716505e7904d24ee3730de25310cb82cd3a1e6e71" # nolint
  ),
  dir = "_targets/user/data_raw/ref_aln"
)

# Download and unzip local fern GenBank database (release 251) from figshare
archive::archive_extract(
  contentid::resolve(
    "hash://sha256/ec689bcf9e97328d5aa200559367f894f4c475a342e33e797a87140c7ca372f0" # nolint
  ),
  dir = "_targets/user/data_raw/restez/sql_db"
)
