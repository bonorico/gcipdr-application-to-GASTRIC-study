


## Install 'JohnsonDistribution' dependency (only available on CRAN archives)

if (!require("JohnsonDistribution"))
{
    url <- "https://cran.r-project.org/src/contrib/Archive/JohnsonDistribution/JohnsonDistribution_0.24.tar.gz"
    pkgFile <- "JohnsonDistribution_0.24.tar.gz"
    download.file(url = url, destfile = pkgFile)
    install.packages(pkgs=pkgFile, type="source", repos=NULL)
    unlink(pkgFile)
    library("JohnsonDistribution")
 }

                                        #

if (!require("gcipdr"))
{
    devtools::install_github("bonorico/gcipdr")
    library(gcipdr)
}


if (!require("st4sas")) {
    devtools::install_github("bonorico/st4sas")
    library(st4sas)
}

                                        #

if (!require("meta")) {
    install.packages("meta")
    library(meta)
}

                        #

if (!require("metafor")) {
    install.packages("metafor")
    library(metafor)
}


                        #

if (!require("coxme")) {
    install.packages("coxme")
    library(coxme)
}

                        #

if (!require("surrosurv")) {
    install.packages("surrosurv")
    library(coxme)
}

                                        #

library(survival)
                                        #


if (!require("frailtypack"))
{
    url <- "https://cran.r-project.org/src/contrib/Archive/frailtypack/frailtypack_3.2.0.1.tar.gz"
    pkgFile <- "frailtypack_3.2.0.1.tar.gz"
    download.file(url = url, destfile = pkgFile)
    untar("frailtypack_3.2.0.1.tar.gz",list=TRUE)  ## check contents
    untar("frailtypack_3.2.0.1.tar.gz", files="frailtypack/data/gastadj.rda")
    load("frailtypack/data/gastadj.rda")
    unlink(pkgFile)
    unlink("frailtypack", TRUE)
## install.packages ("frailtypack", dependencies = T, type = "source", repos = "https://cloud.r-project.org")
## library(frailtypack)
##     data(gastadj)
} else
    data(gastadj)



                      #

if (!require("synthpop")) {
    install.packages("synthpop")
    library(synthpop)
}


                      #

if (!require("xtable")) {
    install.packages("xtable")
    library(xtable)
}


#

if (!require("cowplot")) {
    install.packages("cowplot")
    library(cowplot)
}

#

if (!require("gridExtra")) {
    install.packages("gridExtra")
    library(gridExtra)
}

#

if (!require("gridGraphics")) {
    install.packages("gridGraphics")
    library(gridGraphics)
}


#

if (!require("dplyr")) {
    install.packages("dplyr")
    library(dplyr)
}





source("code/datarebuildflex_v2.R")

source("utils/my_models_utils.R")

source("utils/my_utils.R")

source("utils/vae_generator.R")  ## R 4.02: Error: Installation of TensorFlow not found.

