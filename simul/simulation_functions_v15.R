## simulation_functions_v15.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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



#### simulate survival data with correlated random effects
## same process as in v13 but different seeding. Cannot merge.
## masking function in previous version
## NOTE**: to mantain correlation structure in scenario 2 we just rescale beta. Drawing new beta and frailty changes correlation.
## es: in second_scenario do --> re2 <- MASS::mvrnorm(K, c(0, beta[2]), Sigma); frailty2 <- exp(re2[, 1]); loghr2 <- re2[, 2], then cor(frailty2, loghr2) = Sigma, BUT cor(frailty1, loghr2) or cor(frailty2, loghr1) is NOT Sigma anymore (here typically smaller or wrong sign) which results in estimation artifacts in coxme. RESCALE: loghr2 <- loghr1 + diff(beta)  --> equal results as below ***

generateDataScenario <- function(R = 1, K = 2, n = 1000, gam = c(1.5, 1.05), lam = c(0.1, 0.27),
                                 beta = c(-0.5, -0.005), rdn_sd = c(0.7, 0.15), rho = 0.6,
                                 rateC = c(0.03, 0.05), second_scenario = FALSE, seed = 61 )

{
    n <- ifelse(second_scenario, n/2, n)
    simd <- lapply(1:R,
                     function(i)
                     {
#### list of R independent experiments
                         set.seed(seed + i)
                         cov <- rho*(rdn_sd[1]*rdn_sd[2])
                         Sigma <- rbind(c(rdn_sd[1]^2, cov),
                                        c(cov, rdn_sd[2]^2)
                                        )
                         re <- MASS::mvrnorm(K, c(0, 0), Sigma)
                         frailty <- exp(re[, 1])
                         loghr1 <- beta[1] + re[, 2]
                         if (second_scenario)                       ### non proportional hazard
                             loghr2 <- beta[2] + re[, 2]   

                         s <- do.call(
                             "rbind",
                             mclapply(1:K,
                                      function(j)
                                      {
                                          set.seed(seed + i + j)
                                          s1 <- data.frame(
                                              Weibtimes(n,              ### surv times stratified by center
                                                        lam[1],
                                                        gam[1],
                                                        loghr1[j],
                                                        rateC[1],
                                                        frailty[j] ),
                                              trialID = j,
                                              scen = 1
                                          )
                                          s2 <- NULL
                                          if (second_scenario)
                                              s2 <- data.frame(
                                                  Weibtimes(n,
                                                            lam[2],
                                                            gam[2],
                                                            loghr2[j],
                                                            rateC[2],
                                                            frailty[j] ),
                                                  trialID = j,
                                                  scen = 2
                                              )
                                          return(rbind(s1, s2))
                                      }
                                      )
                         )
                     }
                     
                     )
    return(simd)  # list of R simulated datasets
}



## source("load_dependencies.R")
## source("my_models_utils.R")
## source("simulation_functions_v13.R")

## f <- Surv( time, status )~( 1 + trt|as.factor(trialID) ) + trt

## k  <- 10
## dat <- generateDataScenario(R = 30, K = k, second_scenario = T, beta = c(-1.5, 0.5), rho = 0 )

## out <- BootStat(dat, "random", f)

## out$Bag



### changed 28.1.2021 restructured to allow for inclusion of ranef correlation

                                       
TableReg <- function(regres, colname)
{
    names <- c("coef", "sd", "ranefsd2", "ranefsd1", "revcov")  # changed 28.1.2021 to account for ranef corr
    out <- data.frame( "MC_MEAN" = unlist(regres$MC_MEAN[names]),  "MC_SD" = unlist(regres$MC_SD[names]) )
    colnames(out) <- paste(colname, colnames(out), sep = ".")
    return(out)
}

                                         #
TableRes <- function(simres, moreinfo = NULL)
{
    KM <- rbind( data.frame(simres$KMipd, type = "Orig. IPD"),
                do.call(
                    "rbind",
                    lapply(names(simres$KMcurves), function(m)
                        data.frame(simres$KMcurves[[m]], type = m )
                        )
                ),
                data.frame(simres$KMma, type = "Meta an.")
                )
    REG <- cbind(data.frame(Est = c("coef", "sd", "hrranefsd", "frailty", "rho"),
                            TableReg(simres$Coxipd, "IPD")),   # changed 28.1.2021 to account for additional frailty term
                 do.call(
                     "cbind",
                     lapply(names(simres$Coxreg), function(m)
                         TableReg(lapply(simres$Coxreg[[m]], function(x) x$Bag), m)
                         )
                 ),
                 TableReg(simres$Coxma, "MA")
                 )
    if (!is.null(moreinfo))
    {
        KM <- cbind(moreinfo, KM)
        REG <- cbind(moreinfo, REG)
    }
    rm(simres)   # memory clean-up !!!
    return(list(KM = KM, REG = REG))
}








#################################  TO DELETE ###################

## 05.02.2021: u only allow for rho in well modeled scenario here  .... u must allow always ...

## generateDataScenario <- function(R = 1, K = 2, n = 1000, gam = c(1.5, 1.05), lam = c(0.1, 0.27),
##                                  beta = c(-0.5, -0.005), rdn_sd = c(0.7, 0.15), rho = 0.6,
##                                  rateC = c(0.03, 0.05), second_scenario = FALSE, seed = 61 )

## {
##     n <- ifelse(second_scenario, n/2, n)
##     simd <- mclapply(1:R,
##                      function(i)
##                      {
## #### list of R independent experiments
##                          set.seed(seed + i)
##                          if (!is.null(rho) & !second_scenario)
##                          {
##                              cov <- rho*(rdn_sd[1]*rdn_sd[2])
##                              Sigma <- rbind(c(rdn_sd[1]^2, cov),
##                                             c(cov, rdn_sd[2]^2)
##                                             )
##                              re <- MASS::mvrnorm(K, c(0, beta[1]), Sigma)
##                              frailty <- exp(re[, 1])
##                              loghr1 <- re[, 2]
##                          }
##                          else
##                          {
##                              frailty <- rgmom(K, mean = 1, sd = rdn_sd[1])
##                              loghr1 <- rnorm(K, beta[1], rdn_sd[2])  # random log HR scenario 1 (scenarios share same SD but not the mean)
##                              if (second_scenario)
##                                  loghr2 <- rnorm(K, beta[2], rdn_sd[2])  # random log HR scenario 2
##                          }

##                          s <- do.call(
##                              "rbind",
##                              mclapply(1:K,
##                                       function(j)
##                                       {
##                                           set.seed(seed + i + j)
##                                           s1 <- data.frame(
##                                               Weibtimes(n,
##                                                         lam[1],
##                                                         gam[1],
##                                                         loghr1[j],
##                                                         rateC[1],
##                                                         frailty[j] ),
##                                               trialID = j,
##                                               scen = 1
##                                           )
##                                           s2 <- NULL
##                                           if (second_scenario)
##                                               s2 <- data.frame(
##                                                   Weibtimes(n,
##                                                             lam[2],
##                                                             gam[2],
##                                                             loghr2[j],
##                                                             rateC[2],
##                                                             frailty[j] ),
##                                                   trialID = j,
##                                                   scen = 2
##                                               )
##                                           return(rbind(s1, s2))
##                                       }
##                                       )
##                          )
##                      }
                     
##                      )
##     return(simd)  # list of R simulated datasets
## }
