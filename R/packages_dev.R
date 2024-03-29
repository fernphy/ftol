# Packages not loaded when running code,
# but listed here so renv knows to track them
library(conflicted)
library(clustermq) # for running targets in parallel
library(visNetwork) # for visualizing targets plan
library(languageserver) # for R in VS-Code
library(vscDebugger) # for R in VS-Code
library(jsonlite) # for R in VS-Code
library(docthis) # for documenting functions
library(gittargets) # for tracking targes cache with git
library(BiocManager) # for installing Bioconductor packages
library(here) # used in reports
library(devtools) # used by restrez
library(mirai) # parallelization
library(nanonext) # parallelization
library(crew) # parallelization