panda_profile <- Sys.getenv("R_PROFILE_USER", unset = "")
panda_root <- if (nzchar(panda_profile)) {
  dirname(normalizePath(panda_profile, mustWork = TRUE))
} else {
  getwd()
}
source(file.path(panda_root, "renv", "activate.R"))
rm(panda_profile, panda_root)
