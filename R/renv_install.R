# DO NOT RUN INTERACTIVELY
# only meant for running when building Docker image from Dockerfile
# options(renv.config.pak.enabled = TRUE)
install.packages("renv")
renv::consent(provided = TRUE)
renv::settings$use.cache(FALSE)
renv::init(bare = TRUE)
Sys.setenv(RENV_CONFIG_REPOS_OVERRIDE = "https://cran.r-project.org")
renv::restore()
