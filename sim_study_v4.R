## sim_study_v4.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
## Copyright (C) 2021 Federico Bonofiglio

    ## This Program is free software: you can redistribute it and/or modify
    ## it under the terms of the GNU General Public License as published by
    ## the Free Software Foundation, either version 3 of the License, or
    ## (at your option) any later version.

    ## This Program is distributed in the hope that it will be useful,
    ## but WITHOUT ANY WARRANTY; without even the implied warranty of
    ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ## GNU General Public License for more details.

    ## You should have received a copy of the GNU General Public License
    ## along with This Program.  If not, see <https://www.gnu.org/licenses/>.


### IT EXECUTES CODE TO PRODUCE TABLE 17 SM 

if (!require("devtools"))
{   ### u might need to 'sudo apt-get install libssl-dev libxml2-dev libcurl4-openssl-dev'
    install.packages("devtools")
    library(devtools)
}



if (!require("JohnsonDistribution"))
{
    url <- "https://cran.r-project.org/src/contrib/Archive/JohnsonDistribution/JohnsonDistribution_0.24.tar.gz"
    pkgFile <- "JohnsonDistribution_0.24.tar.gz"
    download.file(url = url, destfile = pkgFile)
    install.packages(pkgs=pkgFile, type="source", repos=NULL)
    unlink(pkgFile)
 }

                                        #

if (!require("gcipdr"))
{
    devtools::install_github("bonorico/gcipdr")
    library(gcipdr)
}


library(st4sas)  #### TODO : create github repository ...

                                        #
if (!require("meta"))
{
    install.packages("meta")
    library(meta)
}

                        #

if (!require("metafor"))
{
    install.packages("metafor")
    library(metafor)
}


                        #

if (!require("coxme")) {
    install.packages("coxme")
    library(coxme)
}

                        #


if (!require("synthpop"))
{
    install.packages("synthpop")
    library(synthpop)
}


                      #

if (!require("xtable")) {
    install.packages("xtable")
    library(xtable)
}


#

if (!require("cowplot"))
{
    install.packages("cowplot")
    library(cowplot)
}

#

if (!require("gridExtra")) {
    install.packages("gridExtra")
    library(gridExtra)
}

#

if (!require("gridGraphics"))
{
    install.packages("gridGraphics")
    library(gridGraphics)
}


library(survival)

library(st4sas)

source("utils/my_models_utils.R")

source("utils/my_utils.R")

source("simul/simulation_functions_v13.R")

source("simul/simulation_functions_v15.R")

### problem with VAE: can't get tf setseed to work in lapply loop--> VAE results will not be reproducible

## source("vae_generator_changed.R")  ## no VAE here

options(mc.cores = 2L)   # to increase parallel cores 

outfolder <- "./results"

### simulation parameters

R <- 30  # fixed

n <- c(50, 500)  # sample size

K <- c(2, 6)  # centers


prms <- list(s1 = list(g = 1.5, l = 0.1, b = -0.5 ),                               # scenario 1
             s2 = list(g = c(1.5, 1.05), l = c(0.1, 0.27), b = c(-0.5, -0.005) ),  # scenario 2
             s3 = list(g = c(1.5, 1.05), l = c(0.1, 0.2), b = c(-1.5, 0.5) ) )     # scenario 3


                                        # cannot use mclapply with keras. It freezes. Run simul separately for vae

## simulation 1: frailty = 0.7 and ranefSD = 0.0

system.time(
    out <- lapply(2:3, function(s)     ### only mis-modeled scenarios
        lapply(1:2, function(i)
        {
            simdR <- generateDataScenario(R,
                            K[i],
                            n[i],
                            prms[[s]]$g,
                            prms[[s]]$l,
                            prms[[s]]$b,
                            rdn_sd = c(0.7, 0),
                            rho = 0,
                            second_scenario = ifelse(length(prms[[s]]$g) > 1, TRUE, FALSE),
                            seed = 20 + s + i
                            )
            ad <- GenerateArtificialDataMultiCenter(simdR, H = 100)

            ad$IPD <- list(simdR)

            return(ad)

        }
            )
            
        )
)


## saveRDS(out, file.path(outfolder, "simulraw_4.rds"))  # in case you save results above

## out <- readRDS(file.path(outfolder, "simulraw_4.rds"))  # in case you save results above

###                data.frame(Scenario = s, "N (Ncenter)" = paste0(n[i], " (",K[i],")") )
#### TEST PH ASSUMPTION ON MIS-MODELED DATA (todo: check log-log km)
system.time(
outtest <- lapply(out, function(s)
    lapply(s, function(m)
        lapply(m[-c(1, 6)], function(r)  ## repetition  c(2, 7)
        {            
            res <- lapply(r, function(d)
                BootStat(d, "residuals", reg_formula = Surv( time, status )~trt + strata(trialID)  )                
                ) 
                                        # return experimental averages 
            return(
                list(Bag = BagInfer(lapply(res, function(x) x$Bag)),
                     BagSE = BagInfer(lapply(res, function(x) x$BagSE)),
                     Bag95QQ = BagInfer(lapply(res, function(x) x$Bag95QQ)) )

            )
        }
        ) ) )
)
                                        #
meths <- c("IPDbootDist", "GCgamma", "GCjohn", "CART", "IPD")

outtab <- do.call("rbind", lapply(1:2, function(s)
    do.call("rbind",
            lapply(1:2, function(i)
                data.frame(Scenario = s,
                           "N (Ncenter)" = paste0(n[i], " (",K[i],")"),
                           do.call("cbind",
                                   lapply(meths,
                                          function(m){                                   
                                              q <- outtest[[s]][[i]][[m]]$Bag$chisq
                                              p <- pchisq(q, 1, lower = F)
                                              o <- data.frame(c(q, p))
                                              rownames(o) <- c("chisq", "pval")
                                              return(o)
                                          }
                                          )
                                   )
                           )
                           )
                )
            )
    )

names(outtab)[-c(1:2)] <- meths 

outtab <- outtab[, c(1:2, 7, 3, 6, 4:5 )]
outtab <- data.frame(outtab[, c(1:2)], rep(c("$\\chi^2$", "$p$"), 4), outtab[, -c(1:2)])
names(outtab)[3] <- "est" 

chis <- apply(round(outtab[which(outtab$est != "$p$"), -c(1:3)], 5), 2, as.character)

outtab_latex <- outtab
outtab_latex[which(outtab_latex$est == "$p$"), -c(1:3)] <- t(apply(outtab[which(outtab$est == "$p$"), -c(1:3)], 1,
                                                       function(x)
                                                           Boldnums(x, TRUE, (x <= 0.05), 5 )
                                                       )
                                                 )

outtab_latex[which(outtab_latex$est != "$p$"), -c(1:3)] <- chis

outtab_latex <- make.multirow(outtab_latex, 1)

outtab_latex[, 2] <- rep(c("\\multirow{2}{4em}{50 (2)}", " ", "\\multirow{2}{4em}{500 (6)}", " "), 2)                                        

outsv <- list(tab_ph_sim = outtab, tab_ph_sim_latex = outtab_latex)

saveRDS(outsv, file.path(outfolder, "torep_simul_4.rds") )

## outsv <- readRDS(file.path(outfolder, "torep_simul_4.rds") )




### is actually nelson-aalen
system.time(
    outKM <- do.call("rbind",
                     lapply(1:2, function(s)
                         do.call("rbind",
                                 lapply(1:2, function(i)
                                     do.call("rbind",
                                             lapply(meths, function(m)
                                                 data.frame(Scenario = s,
                                                            "N (Ncenter)" = paste0(n[i], " (",K[i],")"),
                                                            type = m,
                                                            AverageKM(lapply(out[[s]][[i]][[m]], function(d)
                                                                BootStat(d, "km",
                                                                         Surv( time, status )~trt + strata(trialID),
                                                                         cumhaz = TRUE )$Bag
                                                                )
                                                                )
                                                            )
                                                 )
                                             )
                                     )
                                 )
                         )
                     )
)

outKM$km <- log(outKM$km)  ### log BH
outKM$type <- as.factor(as.character(outKM$type))
   
levels( outKM$type)[-1] <- c("GC few. mom.", "GC all mom.", "Orig. IPD", "Dist. boot.")
outKM$type <-  relevel(outKM$type, "Orig. IPD")
outKM$Scenario <- factor(outKM$Scenario, labels = c("Mis-modeled-1", "Mis-modeled-2") )
outKM$trt <- factor(outKM$trt, labels = c( "Control", "Treatment") )
   
# 

## plot data TODO u can just use KM plots -og-log transformed from v1

GCgam_50_2 <- subset(outKM, N..Ncenter. == "50 (2)" & (type == "Orig. IPD" | type == "GC few. mom.")) 
GCgam_500_6 <- subset(outKM, N..Ncenter. == "500 (6)" & (type == "Orig. IPD" | type == "GC few. mom.")) 

GCjohn_50_2 <- subset(outKM, N..Ncenter. == "50 (2)" & (type == "Orig. IPD" | type == "GC all mom.")) 
GCjohn_500_6 <- subset(outKM, N..Ncenter. == "500 (6)" & (type == "Orig. IPD" | type == "GC all mom.")) 

CART_50_2 <- subset(outKM, N..Ncenter. == "50 (2)" & (type == "Orig. IPD" | type == "CART")) 
CART_500_6 <- subset(outKM, N..Ncenter. == "500 (6)" & (type == "Orig. IPD" | type == "CART")) 

Distboot_50_2 <- subset(outKM, N..Ncenter. == "50 (2)" & (type == "Orig. IPD" | type == "Dist. boot.")) 
Distboot_500_6 <- subset(outKM, N..Ncenter. == "500 (6)" & (type == "Orig. IPD" | type == "Dist. boot."))


 ### plot KM (side by side)

mp <- aes( x= time, y= km, color = trt) 
steps <- geom_step(alpha=0.5, size = 1.3)
panel <- facet_grid( Scenario ~ type)
xl <- xlab("Time");  yl <- ylab("log Nelson-Aalen")
textsize <- theme_bw() + theme(axis.text = element_text(size=18), axis.title = element_text(size=20), strip.text.x = element_text(size = 19), strip.text.y = element_text(size = 19), legend.title= element_blank(), legend.text=element_text(size=18), legend.position="bottom" )



pdf(file.path(outfolder, "bhPHsimul_sidebyside_GCgam_50_2.eps"), width = 10.4, height = 8.3)

ggplot(GCgam_50_2, mp) + steps + panel + xl + yl + textsize 

dev.off()

#
pdf(file.path(outfolder, "bhPHsimul_sidebyside_GCgam_500_6.eps"), width = 10.4, height = 8.3)

ggplot(GCgam_500_6, mp) + steps + panel + xl + yl + textsize 

dev.off()

#

pdf(file.path(outfolder, "bhPHsimul_sidebyside_GCjohn_50_2.eps"), width = 10.4, height = 8.3)

ggplot(GCjohn_50_2, mp) + steps + panel + xl + yl + textsize

dev.off()

#
pdf(file.path(outfolder, "bhPHsimul_sidebyside_GCjohn_500_6.eps"), width = 10.4, height = 8.3)

ggplot(GCjohn_500_6, mp) + steps + panel + xl + yl + textsize

dev.off()

#
pdf(file.path(outfolder, "bhPHsimul_sidebyside_CART_50_2.eps"), width = 10.4, height = 8.3)

ggplot(CART_50_2, mp) + steps + panel + xl + yl + textsize 

dev.off()

#
 pdf(file.path(outfolder, "bhPHsimul_sidebyside_CART_500_6.eps"), width = 10.4, height = 8.3)

ggplot(CART_500_6, mp) + steps + panel + xl + yl + textsize 

dev.off()

                                        
#

pdf(file.path(outfolder, "bhPHsimul_sidebyside_Distboot_50_2.eps"), width = 10.4, height = 8.3)

ggplot(Distboot_50_2, mp) + steps + panel + xl + yl  + textsize 

dev.off()

#
pdf(file.path(outfolder, "bhPHsimul_sidebyside_Distboot_500_6.eps"), width = 10.4, height = 8.3)

ggplot(Distboot_500_6, mp) + steps + panel + xl + yl  + textsize 

dev.off()


