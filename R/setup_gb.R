library(restez)

# Specify location to download GenBank database
restez_path_set("/data_raw")

# Download plant database
db_download(preselection = 1)

# Create database
restez_connect()
db_create(min_length = 10, max_length = 200000)