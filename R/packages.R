# Specify packages used in workflow
workflow_packages <- c(
  "conflicted",
  "assertr",
  "glue",
  "ape",
  "fs",
  "future",
  "future.callr",
  "taxastand",
  "tflow",
  "tidyverse",
  "contentid"
)

# Load the packages for interactive sessions
if (interactive()) invisible(
  lapply(
    workflow_packages, library, character.only = TRUE, quietly = TRUE)
  )
