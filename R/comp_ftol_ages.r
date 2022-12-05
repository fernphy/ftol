# Compare ages of major mononphyletic clades across FTOL versions
#
# Compares between current version in {ftolr} ("ftol") and the one built by
# _targets.R ("new_ftol")

# Load packages ----

library(tidyverse)
library(ftolr)

source("_targets.R")

# Define functions -----

#' Make a tibble of dates and their node numbers in a phylogeny
#'
#' @param phy Phylogeny with dates as a number at each node
#'
#' @return Tibble of dates
make_ages_tbl <- function(phy) {
  ages_values = tibble(bs = phy$node.label) %>%
    pull(bs)

  tibble(
    # hard-code node ID: internal nodes start after tip nodes,
    # and phy$node.label is in the same order as internal nodes
    node = 1:Nnode(phy) + Ntip(phy),
    age = ages_values)
}

#' Relabel a tree with ages at each internal node
#'
#' Should only be used with ultrametric (dated) trees
#'
#' @param tree List of class 'phylo'; the phylogenetic tree
#' @param decimals Numeric vector of length 1; the number of decimals
#'   to round the ages
#' @return Phylogenetic tree with internal nodes labeled by their ages
#' @noRd
label_with_ages <- function(tree, decimals = NULL) {
  # Get total height (age) of tree
  total_height <- max(ape::node.depth.edgelength(tree))
  # Get length for each node from root
  edge_length_all <- ape::node.depth.edgelength(tree)
  # Subset to only internal nodes
  # (the first n edge lengths are tips, the rest are internal nodes)
  edge_length_internal <- edge_length_all[-seq_len(ape::Ntip(tree))]
  # Convert branch lengths to time since root
  edge_age_internal <- total_height - edge_length_internal
  if (!is.null(decimals)) {
    edge_age_internal <- round(edge_age_internal, decimals)
  }
  # Relabel nodes with ages
  tree$node.label <- edge_age_internal
  return(tree)
}

# Load data ----

tar_load(con_monophy_by_clade)
tar_load(sanger_con_tree_dated)

# Set taxonomic levels to compare
taxa_levels_check <- c("family", "suborder", "order")

# Process previous FTOL data -----
ftol <- ft_tree()
ftol_dates <- ft_tree(label_ages = TRUE) %>%
  make_ages_tbl()

ftol_taxonomy <- ftol_taxonomy %>%
  filter(species != "Zygnema_circumcarinatum")

ftol_mono_test <- assess_monophy(
    taxon_sampling = ftol_taxonomy,
    tree = ftol,
    tax_levels = taxa_levels_check
)

ftol_monophy_by_clade = map_df(
    seq_along(taxa_levels_check),
    ~get_result_monophy(ftol_mono_test, .)
  )

ftol_age <-
  ftol_monophy_by_clade %>%
  filter(monophyly == "Yes") %>%
  transmute(taxon, node = as.numeric(mrca)) %>%
  left_join(ftol_dates, by = "node") %>%
  select(taxon, ftol_age = age)

# Process new FTOL data ----
ftol_new <- tar_read(sanger_con_tree_dated)

ftol_new_dates <- ftol_new %>%
  label_with_ages() %>%
  make_ages_tbl()

ftol_new_taxonomy <- tar_read(sanger_sampling) %>%
  filter(species != "Zygnema_circumcarinatum")

ftol_new_mono_test <- assess_monophy(
    taxon_sampling = ftol_new_taxonomy,
    tree = ftol_new,
    tax_levels = taxa_levels_check
)

ftol_new_monophy_by_clade = map_df(
    seq_along(taxa_levels_check),
    ~get_result_monophy(ftol_new_mono_test, .)
  )

new_ftol_age <-
  ftol_new_monophy_by_clade %>%
  filter(monophyly == "Yes") %>%
  transmute(taxon, node = as.numeric(mrca)) %>%
  left_join(ftol_new_dates, by = "node") %>%
  select(taxon, new_ftol_age = age)

# Check difference in ages ----
inner_join(ftol_age, new_ftol_age, by = "taxon") %>%
  mutate(diff = ftol_age - new_ftol_age) %>%
  arrange(diff) %>%
  print(n = 100)
