# RP package template version >= 0.0.77

## Uncomment following code to install package from source directory
# rp::installAndVerify(
#   install = devtools::install,
#   requirement = sprintf("== %s", desc::desc_get_version("/mnt/vol/package_source_dir")),
#   package = "/mnt/vol/package_source_dir"
# )

## Uncomment following code to install package(s) directly from Bitbucket repository
## SSH keys should be copied to container before installation
# ssh_keys <- git2r::cred_ssh_key(file.path("/home/rstudio/.ssh/id_rsa.pub"), file.path("/home/rstudio/.ssh/id_rsa"))
# rp::installAndVerify(
#     install = devtools::install_git,
#     url = "ssh://git@stash.intranet.roche.com:7999/some/repo/url.git",
#     ref = "master", # branch name or commit hash
#     credentials = ssh_keys,
#     package = "package_name",
#     requirement = "*",
#     subdir = NULL, # provide if package is in subdir
#     upgrade = "never"
#   )
