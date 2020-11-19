Fern Tree of Life (FTOL) Data Release v0.0.1 README
===================================================

This README.md file was generated on 18 November, 2020

The goal of the FTOL project is to generate a maximally sampled phylogenetic
tree for all extant fern species. It will be continuously updated as new data
become available. Each release of the tree and associated metadata has a version
number available on the [project github
repository](https://github.com/joelnitta/ftol). The current data release
includes GenBank accessions with allowed dates from 1980-01-01 to 2020-06-30.

Files
-----

Files in this data release include:

-   `ftol_plastid_accs.csv`: Table of GenBank accessions used for phylogenetic
    analysis. Columns include `sci_name` (scientific name including author after
    taxonomic standardization), `species` (same as `sci_name` but without
    author), a series of columns named after genes (e.g., `atpA`, `rbcL`), and
    one column named `full_plastome`. The columns named after genes and
    `full_plastome` each contain the GenBank accession for that gene or full
    plastome. 4,510 rows. MD5 hash: cecd1d95bf1ae3f1a84e429f980d0458

-   `ftol_plastid_concat.fasta`: Concatenated alignment of plastid genes used
    for phylogenetic analysis. 4,510 species x 48,167 base pairs, including 62
    genes. 93.23% missing data. FASTA format. MD5 hash:
    9fa13e54abdfd5c562d1f6a5ee411fe9

-   `ftol_plastid_parts.csv`: Table of gene partitions in concatenated
    alignment. Columns include `gene` (gene name), `start` (gene start position,
    inclusive), `end` (gene end position, inclusive). 62 rows. MD5 hash:
    771718a5177ad34013e7e2a995635f09

-   `ftol_plastid.tre`: Fern tree of life (FTOL) based on plastid genes, not
    dated. 4,510 species (one tip per species), including 4,492 ingroup and 18
    outgroup taxa. Consensus maximum-likelihood tree. Branch lengths in units of
    expected change per site. Values at nodes are support values as measured
    with ultra-fast bootstrap approximation. Newick format. MD5 hash:
    609d4f6f34b8964992d7bd8175c2b0ea

-   `ftol_plastid_dated.tre`: Fern tree of life (FTOL) with divergence times
    estimated on `ftol_plastid.tre` using fossil calibration points. Branch
    lengths in millions of years before present. Newick format. MD5 hash:
    ac2f6da3f1581d09347678eb6e7f28a0

Methods Summary
---------------

The approach leverages the fact that only a handful of genes have been
intensively sequenced in ferns, and account for the vast majority of sequences
on GenBank: *rbcL*, *atpA*, *atpB*, and *rps4*. These four genes (ca. 5 kb
total) provide maximum phylogenetic breadth, but lack the resolving power
necessary to infer a robustly supported tree across all ferns. Therefore, a set
of single-copy genes are also sampled from fern plastomes, and combined with the
Sanger sequences to create a mostly empty supermatrix. The plastome data provide
a well-supported backbone, which is filled in by the much more diverse Sanger
sequences.

The workflow is implemented in R v4.0.0 (R Core Team 2020) with the Drake
package v7.12.2 (Landau 2018). After downloading plastid genes and genomes,
taxonomic name resolution (detection of synonyms) is carried out using the
[Catalog of Life (COL)](https://www.catalogueoflife.org/) as a taxonomic
standard. Next, rogue taxa are detected by an all-by-all BLAST search (Camacho
et al. 2009); any sequences whose top three hits match a different family are
excluded. One specimen per species is selected prioritizing those that come from
the same specimen with maximum combined length. Genes are then aligned
individually using MAFFT (Katoh et al. 2002), trimmed with TrimAl
(Capella-Gutiérrez, Silla-Martínez, and Gabaldón 2009), then concatenated into a
supermatrix. Only sequences originating from the same specimen are concatenated.
The tree is inferred with IQTREE (Nguyen et al. 2015), and dated with treePL
(Smith and O’Meara 2012) using 27 fossil calibration points.

All code is available at
<a href="https://github.com/joelnitta/ftol" class="uri">https://github.com/joelnitta/ftol</a>.

A docker image to run the code is available at
<a href="https://hub.docker.com/r/joelnitta/ftol" class="uri">https://hub.docker.com/r/joelnitta/ftol</a>

References
----------

Camacho, C, G Coulouris, V Avagyan, N Ma, J Papadopoulos, K Bealer, and T
Madden. 2009. “BLAST+: Architecture and Applications.” *BMC Bioinformatics* 10
(1): 421.

Capella-Gutiérrez, Salvador, José M. Silla-Martínez, and Toni Gabaldón. 2009.
“trimAl: A Tool for Automated Alignment Trimming in Large-Scale Phylogenetic
Analyses.” *Bioinformatics* 25 (15): 1972–3.
<https://doi.org/10.1093/bioinformatics/btp348>.

Katoh, Kazutaka, Kazuharu Misawa, Keiichi Kuma, and Takashi Miyata. 2002.
“MAFFT: A Novel Method for Rapid Multiple Sequence Alignment Based on Fast
Fourier Transform.” *Nucleic Acids Research* 30 (14): 3059–66.
<https://doi.org/10.1093/nar/gkf436>.

Landau, William Michael. 2018. “The Drake R Package: A Pipeline Toolkit for
Reproducibility and High-Performance Computing.” *The Journal of Open Source
Software* 3 (21): 550. <https://doi.org/10.21105/joss.00550>.

Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh.
2015. “IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating
Maximum-Likelihood Phylogenies.” *Molecular Biology and Evolution* 32 (1):
268–74. <https://doi.org/10.1093/molbev/msu300>.

R Core Team. 2020. *R: A Language and Environment for Statistical Computing*.
Vienna, Austria: R Foundation for Statistical Computing.
<https://www.R-project.org/>.

Smith, Stephen A, and Brian C. O’Meara. 2012. “treePL: Divergence Time
Estimation Using Penalized Likelihood for Large Phylogenies.” *Bioinformatics*
28 (20): 2689–90. <https://doi.org/10.1093/bioinformatics/bts492>.
