library(gert)
library(gittargets)
library(assertthat)
library(targets)
library(fs)

# Update ftol_data (**output** of ftol workflow), commit,
# and push to github

# Check targets status
targets_status <- tar_git_status_targets()

assert_that(
  nrow(targets_status) == 0,
  msg = "One or more outdated targets detected. Fix, then run tar_git_snapshot()") # nolint

# Check data (_targets store) status
data_status <- tar_git_status_data()

assert_that(
  nrow(data_status) == 0,
  msg = "Data repo not clean. Fix, then run tar_git_snapshot()")

# Check code status
code_status <- tar_git_status_code()

assert_that(
  nrow(code_status) == 0,
  msg = "Code repo not clean. Fix, then run tar_git_snapshot()")

# Write out cc0 license
readr::write_lines(
  cc0(),
  file = "ftol_data/LICENSE")

# Add all files in ftol_data
ftol_data_files <- list.files("ftol_data")

git_add(
  files = ftol_data_files,
  repo = "ftol_data"
)

# Format message:
# commit hash of code repo
# plus comment of code repo
msg <- glue::glue("code={git_commit_info()$id}
({git_commit_info()$message})")

# Commit
git_commit(
  repo = "ftol_data",
  message = msg
)

# Push
