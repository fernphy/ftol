# Register input data for contentid
# Registration only needs to be done once (ever)

library(contentid)
conflicted::conflict_prefer("register", "contentid")

# - Patel et al. 2019 SI
register("https://github.com/susanfawcett/Thelypteridaceae_Spore_Phylogeny/files/8794160/Patel_2019_SupplementaryInfo_1.xlsx") # nolint
# "hash://sha256/233607dc3945dc0f764c44d1171f8bd8bdfe50c4028c9c44e82965e5a5f11fdc" # nolint

# - Testo and Sundue 2016 SI
register("https://github.com/wtesto/4000fernPhylogeny/raw/main/supplementaryTables.xlsx") # nolint
# "hash://sha256/c629f4617e7e2329a10cb1b207b82a6720653a67eafcaaa97cc6ee891ae7fdf7" # nolint

# - NCBI taxonomy v 2022-05-01
register("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2022-05-01.zip") # nolint
# "hash://sha256/883f9d06034178602ada5cff7790903495e9d8e89860aedef7749f931f9c5a23" # nolint

# - fern fossil data: fern_fossils.csv v1.0.0
register("https://raw.githubusercontent.com/fernphy/ferncal/6beb6b59c007d94d58f35658732ff561a9d6a537/fern_fossils.csv") # nolint
# "hash://sha256/55fd2f21d8e26e4604d9128871f9435ede08f75efc8ae64ce56c671f8d605a1e" # nolint

# - reference alignments on FigShare: ref_aln.tar.gz
register("https://figshare.com/ndownloader/files/34604057")
# "hash://sha256/3c37fb9478a8d6d1d5cf12652f01e04c3187db64923be824ca689d924facde18" #nolint