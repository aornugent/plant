#!/bin/bash
set -euo pipefail

# Requires these hosts in the environment's allowed domains:
#   cloud.r-project.org, packagemanager.posit.co, github.com

# Base image ships a deprecated ondrej/php PPA whose metadata changed; drop it
rm -f /etc/apt/sources.list.d/*ondrej*

apt-get update --allow-releaseinfo-change
apt-get install -y --no-install-recommends --fix-missing \
  wget ca-certificates libcurl4-openssl-dev libssl-dev \
  libxml2-dev cmake g++ make libboost-dev gh

# CRAN R 4.5 repo for Noble (requires cloud.r-project.org in allowed domains).
# Not -q: a quiet failure here aborts the script with no output at all.
wget -O /etc/apt/trusted.gpg.d/cran_r45.asc \
  https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc
echo "deb https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" \
  > /etc/apt/sources.list.d/cran-r45.list

apt-get update
apt-get install -y --no-install-recommends --fix-missing r-base r-base-dev

Rscript -e '
options(
  repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"),
  HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
    paste(getRversion(), R.version$platform, R.version$arch, R.version$os)),
  Ncpus = parallel::detectCores()
)
pkgs <- c(
  "arrow", "dplyr", "fs", "nabor", "pkgload",
  "purrr", "qrng", "secretbase", "storr", "tibble",
  "thor", "callr", "crew", "testthat", "remotes",
  "Rcpp", "BH", "R6", "crayon", "magrittr",
  "rlang", "tidyselect", "tidyr", "roxygen2", "uuid",
  "numDeriv", "knitr", "rmarkdown", "ggplot2", "logger",
  "withr", "pkgbuild", "devtools"
)
install.packages(pkgs)
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) {
  message("Failed to install: ", paste(missing, collapse = ", "))
  quit(status = 1)
}
# github.com archive URLs are 403 for repos outside the session; clone instead.
remotes::install_git("https://github.com/richfitz/RcppR6.git", upgrade = "never")
remotes::install_git("https://github.com/traitecoevo/odelia.git", upgrade = "never")
for (p in c("RcppR6", "odelia")) {
  if (!requireNamespace(p, quietly = TRUE)) quit(status = 1)
}
'
