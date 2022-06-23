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

# - fern fossil data: fern_fossils.csv v1.0.1
register("https://raw.githubusercontent.com/fernphy/ferncal/34a16347db3cf068bbe46d3ecfad78ac98b0676e/fern_fossils.csv") # nolint
# "hash://sha256/153139fbb560442ad46770a04be370f5884cd6396e0eaf05d6875513f80d072c" # nolint

# - reference alignments on FigShare: ref_aln.tar.gz
register("https://figshare.com/ndownloader/files/36005603")
# "hash://sha256/388b53201a8626d4b41851e716505e7904d24ee3730de25310cb82cd3a1e6e71" #nolint

# - restez db on FigShare: restez_sql_db.tar.gz
register("https://figshare.com/ndownloader/files/36005606")
#  "hash://sha256/8059a845c6570eeffb6fe08c29e178a9dc223ab6f929a1b6c6b374e160f21410" #nolint