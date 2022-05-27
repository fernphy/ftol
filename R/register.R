# Register local input data for contentid
# Registration only needs to be done once (ever)

library(contentid)
conflicted::conflict_prefer("register", "contentid")

# local registry ----

# Local data files that aren't available via download (or FTP files)
# need to be registered to a local file (./utils/local.tsv).

# - Patel et al. 2019 SI
register("_targets/user/data_raw/44201appS1.xlsx",
  registries = "utils/local.tsv")
# "hash://sha256/5bef560a9c02e1fb99ea15f008ea371de49fef0d29a4528d48ebcca4dad9cfc0" # nolint

# - Testo and Sundue 2016 SI
register("_targets/user/data_raw/1-s2.0-S1055790316302287-mmc2.xlsx",
  registries = "utils/local.tsv")
# "hash://sha256/3438efbd1fbc3513dd6bebe2cd474f59f8db0009bd6cda290190f5b5364ee6b0" # nolint

# online registry -----

# - NCBI taxonomy v 2022-05-01

register("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2022-05-01.zip") # nolint
# "hash://sha256/883f9d06034178602ada5cff7790903495e9d8e89860aedef7749f931f9c5a23" # nolint

# - fern fossil data

# fern_fossils.csv v1.0.0
register("https://raw.githubusercontent.com/fernphy/ferncal/6beb6b59c007d94d58f35658732ff561a9d6a537/fern_fossils.csv") # nolint
# "hash://sha256/55fd2f21d8e26e4604d9128871f9435ede08f75efc8ae64ce56c671f8d605a1e" # nolint

# reference alignments on FigShare

# ref_aln.tar.gz
register("https://figshare.com/ndownloader/files/34604057")
# "hash://sha256/3c37fb9478a8d6d1d5cf12652f01e04c3187db64923be824ca689d924facde18" #nolint