library(contentid)

# Register local input data for contentid
#
# Local data files that aren't available via download need to be registered
# to a local file (./utils/local.tsv).
# Registration only needs to be done once (ever)

register("_targets/user/data_raw/Patel_2019_SI.xlsx",
  registries = "utils/local.tsv")
