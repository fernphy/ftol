# Raw data README

Raw data files used for Fern Tree of Life (FTOL) project

- 1-s2.0-S1055790316302287-mmc2.xlsx: Supplemental information from 
Testo and Sundue 2016 Mol. Phy. Evo. 105: 200-211. 
Downloaded from https://doi.org/10.1016/j.ympev.2016.09.003

- accs_exclude.csv: GenBank accessions to exclude from analysis.

- equisetum_subgenera.csv: Table with Equisetum species and their subgenera. Used
for molecular dating analysis.

- fern_fossils.csv: Ages of fern fossils used for molecular dating analysis.

- plastome_outgroups.csv: Spreadsheet of plastome accessions used as outgroups in
plastid phylogenetic analysis. Columns: group name, species, acccession.

- ppgi_taxonomy_mod.csv: Spreadsheet of taxonomic system for ferns and lycophytes
at genus level and higher following Pteridophyte Phylogeny Group (PPG) I 2016,
modified with new genera added since that time and additional genera used in the
World Ferns list but not included in PPG I.

- ref_aln: Folder containing reference FASTA files for extracting sequences with 
SuperCrunch. Each FASTA file is named after the target locus (gene or spacer
region). Sequences within each FASTA file named by species and GenBank accession
number separated by an underscore. All reference FASTA files produced by
prep_ref_seqs_plan.R

- ref_aln.tar.gz: Same as ref_aln, but compressed and archived.

- target_coding_genes.txt: List of 77 coding genes to extract from plastome
sequences. Based on list from Wei et al 2017 Genome Bio. Evol. 1647-1657,
excluding genes with multiple copies.

- taxdmp_2022-02-01.zip: NCBI taxonomy database. Downloaded from
https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2022-02-01.zip

For a list of hashes generated with md5sum, see hash.txt

## References

Portik, Daniel M., and John J. Wiens. "SuperCRUNCH: A Bioinformatics
toolkit for creating and manipulating supermatrices and other large phylogenetic
datasets". Methods in Ecology and Evolution 11 (2020): 763–72.
https://doi.org/ggx588.

Testo, Weston, and Michael Sundue. "A 4000-species dataset provides new insight
into the evolution of ferns." Molecular Phylogenetics and Evolution 105 (2016):
200-211

Wei, Ran, et al. "Plastid phylogenomics resolve deep relationships among
eupolypod II ferns with rapid radiation and rate heterogeneity." Genome biology
and evolution 9.6 (2017): 1646-1657.
