# FTOL

Code repository to generate an automatically updated **Fern Tree of Life**.

Please see the accompanying paper: Nitta JH, Schuettpelz E, Ramírez-Barahona S, Iwasaki W. An open and continuously updated fern tree of life (FTOL)

All code is in **R**, and workflow is controlled with the [targets](https://github.com/ropensci/targets) package.

## Docker

A [docker image](https://hub.docker.com/repository/docker/joelnitta/ftol) is available to run the code.

## Running the code

To run the analysis, execute `targets::tar_make()` from the root of the repo.
