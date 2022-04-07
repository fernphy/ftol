library(contentid)

# register FTOL input data on figshare
#
# obtain URLs by navigating to
# https://doi.org/10.6084/m9.figshare.19474316 #nolint
# in broswer, and copying URL from "Download" link
# for each file
#
# make sure to select a different version (DOI) if needed
#
# registration only needs to be done once (ever)

# DOI: 10.6084/m9.figshare.19474316.v1 ----

# accs_exclude.csv
register("https://figshare.com/ndownloader/files/34604039")
# "hash://sha256/e8e215870706a03e74f21a7e117246d77ff1029d91599cda53ec14ea7fbcc1ab" #nolint

# equisetum_subgenera.csv
register("https://figshare.com/ndownloader/files/34604045")
# "hash://sha256/a93ec0663a65d687921af3c412279034786fba769d73408c432bd9b738bd37ad" #nolint

# plastome_outgroups.csv
register("https://figshare.com/ndownloader/files/34604051")
# "hash://sha256/36bf35dbe61c4f133ba5b7112681316b4338480fb6b1299b762f893d5e89c6d1" #nolint

# ppgi_taxonomy_mod.csv
register("https://figshare.com/ndownloader/files/34604054")
# "hash://sha256/d94cf3b3230a4fafaadf76b355a9d989cc1645467aab47934a73cba2920fff3f" #nolint

# target_coding_genes.txt
register("https://figshare.com/ndownloader/files/34604060")
# "hash://sha256/304cd16b67e1d4f180624bed3b683c848c42b6ad5d8250cda1f0425e58831ccf" #nolint

# ref_aln.tar.gz
register("https://figshare.com/ndownloader/files/34604057")
# "hash://sha256/3c37fb9478a8d6d1d5cf12652f01e04c3187db64923be824ca689d924facde18" #nolint