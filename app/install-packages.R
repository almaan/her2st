#!/usr/bin/Rscript
require(yaml)

## Installation Script for GeneTech course

# get installed packages
i.pkgs <- rownames(installed.packages())

# packages on CRAN
a.pkgs <- read_yaml("install.yaml")
n.pkgs <- a.pkgs[["CRAN"]]
b.pkgs <- a.pkgs[["BioConductor"]]

# to collect stats
success <- c()
fail <- c()

# install CRAN based packages
if (!(is.null(n.pkgs))) {
  for (p in n.pkgs) {
    # check if already installed
    if (!(p %in% i.pkgs)){
      # try-catch install
      stat <- try(install.packages(p,
                                  dependencies = TRUE,
                                  repos = "https://ftp.acc.umu.se/mirror/CRAN/"))

      # add status
      if (!(is.null(attr(stat,"class")))) {
        fail <- c(fail,p)
        } else {
          success <- c(success,p)
        }
    } else {
      success <- c(success,p)
    }
  }
}

# install Bioconductor packages
if (!is.null(b.pkgs)) {
  for (p in b.pkgs) {
    # check if already installed
    if (!(p %in% i.pkgs)){
      # try-catch install
      stat <- try(BiocManager::install(p))
      # add status
      if (!(is.null(attr(stat,"class")))) {
        fail <- c(fail,p)
      } else {
        success <- c(success,p)
      }
    } else {
      success <- c(success,p)
    }
  }
}

# log-message construction
success <- ifelse(is.null(success),
                  "NONE\n\n",
                  paste(success,collapse = "\n"))

fail <- ifelse(is.null(fail),
               "NONE\n\n",
               paste(fail,collapse = "\n")
               )

txt <- paste(c("\n----\n",
               "Successful installs:\n",
               "----\n",
               success,
               "\n\n----\n",
               "Failed installs:\n",
               "----\n",
               fail),
               collapse  = "")

# print summary message
cat(txt)
