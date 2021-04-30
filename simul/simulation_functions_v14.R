## simulation_functions_v14.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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


### adaptations to run code in sim_study_v2.R. Some functions from _v13 are now MASKED


#### changed 23.6.20 restructured to allow for inclusion of additional frailty effect in tables

Simulations_stat <- function( x, formula, const.haz = FALSE,
                             frailty = FALSE, KM = FALSE, cumhaz = FALSE )
{
    x <- as.data.frame(x)
    if (!KM)
    {
        if (!frailty)
        {
            if (!const.haz)
                out <- try(
                    coxph( formula , data = x ),
                    silent = TRUE
                )
            else
                out <- try(
                    glm( formula, "poisson", data = x ),
                    silent = TRUE
                )
            if (class(out)[1] == "try-error")
                return( list( coef= NA, sd = NA, vcov = NA,
                             ranefsd2= NA, ranefsd1= NA, ranefsdSE = NA,
                             revcov = NA, ranef = NA, lr = NA, df = NA )
                       )
            covmat <- vcov(out)
            return(
                list( coef= coef(out), sd = sqrt(diag(covmat)),
                     vcov = covmat[lower.tri(covmat)], ranefsd2= NA,
                     ranefsd1= NA, ranefsdSE = NA,
                     revcov = NA, ranef = NA, lr = NA, df = NA
                     )
            )
        }
        else
        {
            out <- try( coxme( formula , data = x ) , silent = TRUE  )
            if (class(out) == "try-error")
                return( list( coef= NA, sd = NA, vcov = NA,
                             ranefsd2= NA, ranefsd1= NA, ranefsdSE = NA,
                             revcov = NA, ranef = NA, lr = NA, df = NA)
                       )
            covmat <- vcov(out)
            recovmat <-  as.matrix(VarCorr(out)[[1]])
            ranefsd <-  sqrt(diag(recovmat))
            return(
                list(coef= fixef(out),
                     sd = sqrt(diag(covmat)),
                     vcov = covmat[lower.tri(covmat)],
                     ranefsd2= ifelse(length(ranefsd) > 1,
                                      ranefsd[2], NA),
                     ranefsd1= ranefsd[1],
                     ranefsdSE = NA,
                     revcov = recovmat[lower.tri(recovmat)],
                     ranef = ranef(out)[[1]],
                     lr = NA, df = NA
                     )
            )
        }
    }
    else
    {
        estim <- ifelse(cumhaz, "cumhaz", "surv")
        out <- try( survfit( formula, x ), silent = TRUE)                 
        if (class(out) == "try-error")
            return(data.frame(km = NA, time = NA,
                              km_l = NA, km_u = NA, trt = c(0,1))
                   )
        if (!is.null(out$strata))
        {
            ncontrol <- out$strata[1]
            ntrt <- out$strata[2]
            ntreat <- (ncontrol + 1):(ncontrol + ntrt)
                                        #  to get cumhaz CIs apply -log(surv) on km_l, km_u
            return(
                rbind(
                    data.frame(
                        km = out[[estim]][1:ncontrol], time = out$time[1:ncontrol],
                        km_l = out$lower[1:ncontrol], km_u = out$upper[1:ncontrol], trt = 0 ),
                    data.frame(
                        km = out[[estim]][ntreat], time = out$time[ntreat],
                        km_l = out$lower[ntreat], km_u = out$upper[ntreat], trt = 1
                    )
                )
            )
        }
        else
        {
            return(data.frame(km = out[[estim]], time = out$time,
                              km_l = out$lower, km_u = out$upper, trt = NA )
                   )
        }
    }

}


                                       

#### changed 23.6.20 restructured to allow for inclusion of additional frailty effect in tables

statma <- function(TE, seTE, hubsnames, effect = "HR",
                   statist = c("random", "fixed"))
{
    statist <- match.arg(statist)
    drop <- which(seTE >= 1e+02)
    if (length(drop) > 0)
        out <- try(
            metagen(TE[-drop], seTE[-drop], hubsnames[-drop],
                    sm = effect, method.tau = "SJ"),
            silent = TRUE)
    else
        out <- try(
            metagen(TE, seTE, hubsnames, sm = effect,
                    method.tau = "SJ"),
            silent = TRUE)
    if (class(out)[[1]] == "try-error")
        return(data.frame(coef= NA, sd = NA, vcov = NA, ranefsd2= NA,
                          ranefsd1= NA, ranefsdSE = NA,
                          revcov = NA, ranef = NA, lr = NA, df = NA )
               )
    else
        return( data.frame( coef= out[[paste0("TE.",statist)]],
                           sd = out[[paste0("seTE.",statist)]],
                           ranefsd2 = ifelse(statist[1] == "random",
                                             out$tau, NA),
                           ranefsd1= NA,
                           ranefsdSE = ifelse(statist[1] == "random",
                                              out$se.tau2, NA),
                           revcov = NA, ranef = NA, lr = NA, df = NA )
               )
}


### changed 23.6.20 restructured to allow for inclusion of additional frailty effect in tables

                                       
TableReg <- function(regres, colname)
{
    names <- c("coef", "sd", "ranefsd2", "ranefsd1")  # changed 23.6.20 to account for additional frailty term
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
    REG <- cbind(data.frame(Est = c("coef", "sd", "hrranefsd", "frailty"),
                            TableReg(simres$Coxipd, "IPD")),   # changed 23.6.20 to account for additional frailty term
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








