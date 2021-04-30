## sim_study_v3.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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


### IT EXECUTES CODE TO PRODUCE TABLE 16 SM 

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

if (!require("coxme"))
{
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

if (!require("xtable"))
{
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

if (!require("gridExtra"))
{
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

source("simul/simulation_functions_v14.R")  ## masking some functions from previous source ..is OK

source("simul/simulation_functions_v15.R")  ## masking some functions from previous source ..is OK

### problem with VAE: can't get tf setseed to work in lapply loop--> VAE results will not be reproducible

source("utils/vae_generator_changed.R")

## options(mc.cores = 2L)   # to increase parallel cores 

outfolder <- "./resultse"

### simulation parameters

R <-  30  # fixed

n <- c(50, 500)  # sample size

K <- c(2, 10)  # centers


prms <- list(s1 = list(g = 1.5, l = 0.1, b = -0.5 ),                               # scenario 1
             s2 = list(g = c(1.5, 1.05), l = c(0.1, 0.27), b = c(-0.5, -0.005) ),  # scenario 2
             s3 = list(g = c(1.5, 1.05), l = c(0.1, 0.2), b = c(-1.5, 0.5) ) )     # scenario 3


                                        # cannot use mclapply with keras. It freezes. Run simul separately for vae

## simulation 3: frailty = 0.7, ranefSD = 0.15, and rho = 0.6

system.time(
    out <- lapply(1:3,
                  function(s)
                      lapply(1:2,
                             function(i)
                                 TableRes(
                                     Simulations(R,
                                                 K[i],
                                                 n[i],
                                                 prms[[s]]$g,
                                                 prms[[s]]$l,
                                                 prms[[s]]$b,
                                                 rdn_sd = c(0.7, 0.15),
                                                 rho = 0.6,
                                                 reg_formula = Surv( time, status )~(1 + trt|as.factor(trialID) ) + trt ,       # no trial stratification needed according to ground-truth data
                                                 second_scenario = ifelse(length(prms[[s]]$g) > 1, TRUE, FALSE),
                                                 seed = 20 + s + i,
                                                 DoVae = FALSE
                                                 ),
                                     data.frame(Scenario = s,
                                                "N (Ncenter)" = paste0(n[i],
                                                                       " (",K[i],")")
                                                )
                                 )
                             )
                  )
)


saveRDS(out, file.path( outfolder, "simul_3.rds"))  # in case you save results above

#  out <- readRDS(file.path( outfolder,"simul_3.rds"))       # in case you load already saved results 

# time plan : 21 h approx on 2 cores



outReg <- do.call("rbind",
                  lapply(out,
                         function(s)
                             do.call(
                                 "rbind",
                                 lapply(s,
                                        function(i)
                                            i$REG)
                             )
                         )
                  )


outKM <- do.call("rbind",
                 lapply(out,
                        function(s)
                            do.call("rbind",
                                    lapply(s,
                                           function(i)
                                               i$KM)
                                    )
                        )
                 )

outKM$type <- as.factor(as.character(outKM$type))

    
## run vae 21 h total block  on 2 cores

system.time(
    outVAE <- lapply(1:3,
                     function(s)
                         lapply(1:2,
                                function(i)
                                    TableRes(
                                        Simulations(R,
                                                    K[i],
                                                    n[i],
                                                    prms[[s]]$g,
                                                    prms[[s]]$l,
                                                    prms[[s]]$b,
                                                    rdn_sd = c(0.7, 0.15),
                                                    rho = 0.6,
                                                    reg_formula = Surv( time, status )~(1 + trt|as.factor(trialID) ) + trt ,       # no trial stratification needed according to ground-truth data
                                                    second_scenario = ifelse(length(prms[[s]]$g) > 1, TRUE, FALSE),
                                                    seed = 20 + s + i,
                                                    DoVae = TRUE,
                                                    DoVaeOnly = TRUE,
                                                    noRefan = TRUE
                                                    ),
                                        data.frame(Scenario = s, "N (Ncenter)" = paste0(n[i], " (",K[i],")") ) ) )  )
  )


saveRDS(outVAE, file.path( outfolder, "simulVAE_3.rds"))   # in case you save results above

#  outVAE <- readRDS(file.path( outfolder, "simulVAE_3.rds")) # in case you load already saved results 


outVAEReg <- do.call("rbind", lapply(outVAE, function(s) do.call("rbind", lapply(s, function(i) i$REG)  ))  )

outVAEKM <- do.call("rbind", lapply(outVAE, function(s) do.call("rbind", lapply(s, function(i) i$KM)  ))  )
outVAEKM <- outVAEKM[which(outVAEKM$type == "VAE"), ]
outVAEKM$type <- as.factor(as.character(outVAEKM$type))

## NOTE TODO : VAE results are pretty bad... maybe try simulations with less midlle nodes 6 --> 4 or 5 ? 

# 
outReg[, 16:17] <- outVAEReg[, 6:7]

outReg <- make.multirow(outReg, c(1))

outReg[, 2] <-  rep( c("\\multirow{5}{4em}{50 (2)}", rep("", 4), "\\multirow{5}{4em}{500 (6)}", rep("", 4)), 3)

simulKM <- rbind(outKM, outVAEKM)



 # make distinction between orig. IPD estimate and all others 

levels( simulKM$type)[2:5] <- c("GC few. mom.", "GC all mom.", "IPD boot.", "Dist. boot.")
simulKM$type <-  relevel(simulKM$type, "Orig. IPD")
simulKM$Scenario <- factor(simulKM$Scenario, labels = c("Well-modeled", "Mis-modeled-1", "Mis-modeled-2") )
simulKM$trt <- factor(simulKM$trt, labels = c( "Control", "Treatment") )
   

##### Approximate "true" time-averaged HR (main manuscript Appendix A) for mis-modeled secanrios, using same simulation parameters as above, except n which goes to inf

ref <- lapply(2:3,
              function(s)
                  unlist(
                      lapply(1:2,
                             function(i)
                             {
                                 ipdR <-  generateDataScenario(R,
                                                               K[i],
                                                               50000,
                                                               prms[[s]]$g,
                                                               prms[[s]]$l,
                                                               prms[[s]]$b,
                                                               rdn_sd = c(0.7, 0.15),
                                                               rho = 0.6,        ## new setting !!!
                                                               second_scenario = ifelse(
                                                                   length(prms[[s]]$g) > 1,
                                                                   TRUE,
                                                                   FALSE),
                                                               seed = 20 + s + i
                                                               )
                                 mean(
                                     unlist(
                                         lapply(ipdR,
                                                function(x)
                                                    time_aver_HR(x,
                                                                 prms[[s]]$b[1],
                                                                 prms[[s]]$b[2]
                                                                 )
                                                )
                                     ),
                                     na.rm = T
                                 )
                             }
                             )
                  )
              )


### add reference vector 

Ref <- c(
    c( prms[[1]]$b, outReg[2, 4], 0.15, 0.7, 0.6),  
    c( prms[[1]]$b, outReg[7, 4], 0.15, 0.7, 0.6),
    c( ref[[1]][1], outReg[12, 4], 0.15, 0.7, 0.6),
    c( ref[[1]][2], outReg[17, 4], 0.15, 0.7, 0.6),
    c( ref[[2]][1], outReg[22, 4], 0.15, 0.7, 0.6),
    c( ref[[2]][2], outReg[27, 4], 0.15, 0.7, 0.6)
)



## boldface results 

round_up <- 2
tab_thres <- 0.05

M <- dim(outReg)[2]-3

outReg_bold_coef <- as.data.frame(
    t(
        apply(
            cbind(Ref,
                  outReg[, -c(1,2,3)][, (1:M)%%2 == 1] ),   ### select only mean
            1,  
            function(x)
                Boldnums(x,
                         tab_thres,
                         abs(x - x[1]),
                         round_up
                         )
        )
    )
) 

outREG <- outReg    ### copy data to edit

outREG[, c(F,F,F, (1:M)%%2 == 0)] <- round(outReg[, c(F,F,F, (1:M)%%2 == 0)], round_up)   ### select only MSE

outREG[, c(F,F,F, (1:M)%%2 == 1)] <- outReg_bold_coef[, -1] 

outREG$Est <- factor(outReg$Est, labels = c("$\\beta$", "$\\xi$", "$\\tau$", "$\\rho$", "SD"))

pair <- matrix( c((1:M)[(1:M)%%2 == 1], (1:M)[(1:M)%%2 == 0]), ncol = 2)
 
outReg_comp <- do.call(
    "cbind",
    lapply(1:(M/2),
           function(i)
               apply(outREG[, -c(1:3)][, pair[i, ]], 1,
                     function(x)
                         paste0(x[1], "(", x[2],")")
                     )
           )
)

colnames(outReg_comp) <- names(outReg[, -c(1:3)])[(1:M)%%2 == 1]

Ref[c(2, 7, 12, 17, 22, 27 )] <- NA
 
outREG <- cbind(outREG[, 1:3], as.character(round(Ref, 2)), outReg_comp)  # fourth column is ref value

outREG[, 2] <- as.factor(  outREG[, 2])
levels(outREG[, 2])[-1] <- paste0("\\parbox[t]{2mm}{\\multirow{5}{*}{\\rotatebox[origin=c]{90}{",paste0(n, " (", K,")"), "}}}")


for (i in names(outREG)[-c(1:4)]) {     ### rename NA as empty space
    outREG[[i]] <- as.factor(outREG[[i]])
    lev <- levels(outREG[[i]])
    lev[which(lev == "NA(NA)")] <- ""
    levels(outREG[[i]]) <- lev
}
      


saveRDS(outREG, file.path(outfolder, "torep_simul_3.rds") )





### inspect closing FUP times

tau <- lapply(2:3,
              function(s)
                  lapply(1:2, function(i)
                  {
                      ipdR <-  generateDataScenario(R, K[i],
                                                    50000, prms[[s]]$g, prms[[s]]$l, prms[[s]]$b,
                                                    rdn_sd = c(0.4, 0.25),
                                                    second_scenario = ifelse(
                                                        length(prms[[s]]$g) > 1,
                                                        TRUE,
                                                        FALSE),
                                                    seed = 20 + s + i )

    
                      apply( do.call(
                          "rbind",
                          lapply(ipdR,
                                 function(x)
                                     data.frame(t1 = max(x$time[x$scen == 1], na.rm = T),
                                                t2 = max(x$time[x$scen == 2], na.rm = T) ) ) ), 2, mean )
                  }
                  )
              )   


 ## NOTE on tau: in scenario 2 t2 (first time span) is a bit longer than t1, and so the opposite in scenario 3



### DISCUSSION (10.02.2021): treatment effect in KM plots is no longer visible under IPD simulations with tau and rho .. 

## plot data DO WE NEED IT or shall we just use data from simul 1 ? TODO 27.01.2021 

GCgam_50_2 <- subset(simulKM, N..Ncenter. == "50 (2)" & (type == "Orig. IPD" | type == "GC few. mom.")) 
GCgam_500_6 <- subset(simulKM, N..Ncenter. == "500 (10)" & (type == "Orig. IPD" | type == "GC few. mom.")) 

GCjohn_50_2 <- subset(simulKM, N..Ncenter. == "50 (2)" & (type == "Orig. IPD" | type == "GC all mom.")) 
GCjohn_500_6 <- subset(simulKM, N..Ncenter. == "500 (10)" & (type == "Orig. IPD" | type == "GC all mom.")) 

CART_50_2 <- subset(simulKM, N..Ncenter. == "50 (2)" & (type == "Orig. IPD" | type == "CART")) 
CART_500_6 <- subset(simulKM, N..Ncenter. == "500 (10)" & (type == "Orig. IPD" | type == "CART")) 

VAE_50_2 <- subset(simulKM, N..Ncenter. == "50 (2)" & (type == "Orig. IPD" | type == "VAE")) 
VAE_500_6 <- subset(simulKM, N..Ncenter. == "500 (10)" & (type == "Orig. IPD" | type == "VAE")) 

MA_50_2 <- subset(simulKM, N..Ncenter. == "50 (2)" & (type == "Orig. IPD" | type == "Meta an.")) 
MA_500_6 <- subset(simulKM, N..Ncenter. == "500 (10)" & (type == "Orig. IPD" | type == "Meta an."))

Distboot_50_2 <- subset(simulKM, N..Ncenter. == "50 (2)" & (type == "Orig. IPD" | type == "Dist. boot.")) 
Distboot_500_6 <- subset(simulKM, N..Ncenter. == "500 (10)" & (type == "Orig. IPD" | type == "Dist. boot."))


 ### plot KM (side by side)

mp <- aes( x= time, y= km, color = trt) 
steps <- geom_step(alpha=0.5, size = 1.3)
panel <- facet_grid( Scenario ~ type)
xl <- xlab("Time");  yl <- ylab("Kaplan-Meier")
xlims <- xlim(0, 2)
textsize <- theme_bw() + theme(axis.text = element_text(size=18), axis.title = element_text(size=20), strip.text.x = element_text(size = 19), strip.text.y = element_text(size = 19), legend.title= element_blank(), legend.text=element_text(size=18), legend.position="bottom" )



 pdf("kmSim_sidebyside_GCgam_50_2_simul_3.eps", width = 10.4, height = 8.3)

   ggplot(GCgam_50_2, mp) + steps + panel + xl + yl + xlims + textsize + ylim(0.75, 1)

 dev.off()

#
 pdf("kmSim_sidebyside_GCgam_500_6_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(GCgam_500_6, mp) + steps + panel + xl + yl + xlims + textsize + ylim(0.75, 1)

 dev.off()

#

pdf("kmSim_sidebyside_GCjohn_50_2_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(GCjohn_50_2, mp) + steps + panel + xl + yl + xlims + textsize + ylim(0.75, 1)

 dev.off()

#
 pdf("kmSim_sidebyside_GCjohn_500_6_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(GCjohn_500_6, mp) + steps + panel + xl + yl + xlims + textsize + ylim(0.75, 1)

 dev.off()

#
pdf("kmSim_sidebyside_CART_50_2_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(CART_50_2, mp) + steps + panel + xl + yl + xlims + textsize + ylim(0.75, 1)

 dev.off()

#
 pdf("kmSim_sidebyside_CART_500_6_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(CART_500_6, mp) + steps + panel + xl + yl + xlims + textsize + ylim(0.75, 1)

 dev.off()

                                        #
pdf("kmSim_sidebyside_VAE_50_2_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(VAE_50_2, mp) + steps + panel + xl + yl + xlims + textsize + ylim(0.75, 1)

 dev.off()

#
 pdf("kmSim_sidebyside_VAE_500_6_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(VAE_500_6, mp) + steps + panel + xl + yl + xlims + textsize + ylim(0.75, 1)

 dev.off()
                                        #

pdf("kmSim_sidebyside_MA_50_2_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(MA_50_2, mp) + steps + panel + xl + yl + xlims + textsize + ylim(0.75, 1)

 dev.off()

#
 pdf("kmSim_sidebyside_MA_500_6_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(MA_500_6, mp) + steps + panel + xl + yl + xlims + textsize + ylim(0.75, 1)

 dev.off()
                                        #

pdf(file.path(outfolder, "kmSim_sidebyside_Distboot_50_2.eps"), width = 10.4, height = 8.3)

ggplot(Distboot_50_2, mp) + steps + panel + xl + yl + xlims + textsize + ylim(0.75, 1)

dev.off()

#
pdf(file.path(outfolder, "kmSim_sidebyside_Distboot_500_6.eps"), width = 10.4, height = 8.3)

ggplot(Distboot_500_6, mp) + steps + panel + xl + yl + xlims + textsize + ylim(0.75, 1)

dev.off()

#


 ### plot KM (overlap)

 mp <- aes( x= time, y= km, color = type)
panel <- facet_grid(Scenario ~ trt)
xlims <- xlim(0, 30)
textsize <- theme_bw() + theme(axis.text = element_text(size=18), axis.title = element_text(size=20), strip.text.x = element_text(size = 19), strip.text.y = element_text(size = 19), legend.title=element_blank(), legend.text=element_text(size=18), legend.position="bottom" )


# png("km_sim_GCgamma_50_2.png", width = 1000, height = 800)
pdf("kmSim_overlap_GCgam_50_2_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(GCgam_50_2, mp) + steps + panel + xl + yl + xlims + textsize

  dev.off()

#

 pdf("kmSim_overlap_GCgam_500_6_simul_3.eps", width = 10.4, height = 8.3)

 ggplot(GCgam_500_6, mp) + steps + panel + xl + yl + xlims + textsize

 dev.off()

                                        #


 pdf("kmSim_overlap_GCjohn_50_2_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(GCjohn_50_2, mp) + steps + panel + xl + yl + xlims + textsize

  dev.off()

#


 pdf("kmSim_overlap_GCjohn_500_6_simul_3.eps", width = 10.4, height = 8.3)

 ggplot(GCjohn_500_6, mp) + steps + panel + xl + yl + xlims + textsize

 dev.off()


                                        #


 pdf("kmSim_overlap_CART_50_2_simul_3.eps", width = 10.4, height = 8.3)

 ggplot(CART_50_2, mp) + steps + panel + xl + yl + xlims + textsize

  dev.off()

#


 pdf("kmSim_overlap_CART_500_6_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(CART_500_6, mp) + steps + panel + xl + yl + xlims + textsize

 dev.off()

 #


 pdf("kmSim_overlap_VAE_50_2_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(VAE_50_2, mp) + steps + panel + xl + yl + xlims + textsize

  dev.off()

#


 pdf("kmSim_overlap_VAE_500_6_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(VAE_500_6, mp) + steps + panel + xl + yl + xlims + textsize

 dev.off()

 #


 pdf("kmSim_overlap_MA_50_2_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(MA_50_2, mp) + steps + panel + xl + yl + xlims + textsize

  dev.off()

#


 pdf("kmSim_overlap_MA_500_6_simul_3.eps", width = 10.4, height = 8.3)

  ggplot(MA_500_6, mp) + steps + panel + xl + yl + xlims + textsize

dev.off()

                                        #
pdf(file.path(outfolder, "kmSim_overlap_Distboot_50_2.eps"), width = 10.4, height = 8.3)

ggplot(Distboot_50_2, mp) + steps + panel + xl + yl + xlims + textsize

dev.off()

#


pdf(file.path(outfolder, "kmSim_overlap_Distboot_500_6.eps"), width = 10.4, height = 8.3)

ggplot(Distboot_500_6, mp) + steps + panel + xl + yl + xlims + textsize

dev.off()
