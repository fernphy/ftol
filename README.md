## FTOL

The goal of this project is to generate an automatically updated **Fern Tree of Life**.

All code is in R, and workflow is controlled with the [drake](https://docs.ropensci.org/drake/) package.

A [docker image](https://hub.docker.com/repository/docker/joelnitta/ftol) is available to run code.

To run the analysis, execute `drake::r_make()` from the root of the repo.

For more info, [see the wiki](https://github.com/joelnitta/ftol/wiki).
