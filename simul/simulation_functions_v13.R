## simulation_functions_v13.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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


### SIMULATION MODULES (for sim_study.R and sim_study_v2.R)

### survival times from a Weibull dist with frailty term Z

Aweib <- function(t, lambda, rho, beta) exp(beta)*lambda*t^rho
                                        #
Sweib <- function(t, lambda, rho, beta) exp(-Aweib(t, lambda, rho, beta))
                                        #
invWeib <- function(u, lambda, rho, beta) (- log(u) / (lambda * exp(beta) ))^(1 / rho)    
                                        #

Weibtimes <- function(N, lambda, rho, beta, rateC, Z = 1)
{
  # covariate 
  x <- rbinom(N, 1, 0.5)
    
  # Weibull times
  U <- runif(n=N)
  Tstar <- (- log(U) / (lambda * exp(x * beta) * Z ))^(1 / rho)

  # censoring times
  C <- rexp(n=N, rate=rateC)

  # follow-up times and event indicators
  time <- pmin(Tstar, C)
  status <- as.numeric(Tstar <= C)
                                        # data set
    data.frame( time=time, status=status, trt=x)
}

#### DATA SCENARIOs (Weibull distributed surv data)

                                        # R simulation repetitions
                                        # K centers
                                        # n sample size within center
                                        # lam (positive Weib param) the smaller the bigger mean time
                                        # rateC censoring rate typically less than lam

                                        # changed 4.6.20 allow for beta RE. rdn_sd[1] intercept, rdn_sd[2] slope

## NOTE**: beware, drawing new beta in scenario 2 could generate artifact in correlation (e.g. accidentally bigger/smaller than zero here). Here correlation is not focus. Otherwise consider just rescaling beta ..... 12.02.2021 REPLACE THIS FUNC with that in v_15

generateDataScenario <- function(R = 1, K = 2, n = 1000, gam = c(1.5, 1.05), lam = c(0.1, 0.27),
                                 beta = c(-0.5, -0.005), rdn_sd = c(0.7, 0.0), rho = NULL, 
                                 rateC = c(0.03, 0.05), second_scenario = FALSE, seed = 61 )
{
    n <- ifelse(second_scenario, n/2, n)
    simd <- mclapply(1:R, function(i){   #### list of R independent experiments
        set.seed(seed + i)
        frailty <- rgmom(K, mean = 1, sd = rdn_sd[1]) # frailty distr  # changed 4.6.20 let it be Normal # exp( rnorm(K, sd = rdn_sd[1]) )
        loghr1 <- rnorm(K, beta[1], rdn_sd[2])  # random log HR scenario 1 (scenarios share same SD but not the mean)
        if (second_scenario)
            loghr2 <- rnorm(K, beta[2], rdn_sd[2])  # random log HR scenario 2 (see NOTE**)
        s <- do.call("rbind", mclapply(1:K, function(j)
        {
            set.seed(seed + i + j)
            s1 <- data.frame( Weibtimes(n, lam[1], gam[1], loghr1[j], rateC[1], frailty[j] ), trialID = j, scen = 1 )
            s2 <- NULL
            if (second_scenario)
                s2 <- data.frame( Weibtimes(n, lam[2], gam[2], loghr2[j], rateC[2], frailty[j] ), trialID = j, scen = 2 )
            return(rbind(s1, s2))
        }
        )
        )
    }
    )
    return(simd)  # list of R simulated datasets       
 }



                                        #

GCboot <- function(simd, H, SI_k = 20000, meth = 3, seedAD = 78)
{
    set.seed( seedAD, "L'Ecuyer")     
    out <- try(
        Simulate.many.datasets(list(simd),
                               H, meth, stochastic.integration = TRUE,
                               SI_k = SI_k ),
        silent = TRUE)
    if (class(out[[1]]) == "try-error")
        return( lapply(1:H, function(x) replace(simd[1,], values = NA) )  )
    return(out[[1]][["similar.data"]])   # directly extract artif data

}
                                        #

cart <- function(simd, H, seedAD = 78)
{
    out <- try( syn(simd, seed =  seedAD, m = H ), silent = TRUE )
    if (class(out) == "try-error")
        return( lapply(1:H, function(x) replace(simd[1,], values = NA) )  )
    return(out$syn)

}

                                        #
ChangeVAEarchit <- function(file = file.path(getwd(), "vae_generator.R"),
                            new_visible_dim = 3L, new_intermediate_dim = 6L,
                            source.file = FALSE)
{
    tx <- readLines(file)
    tx2 <- gsub("original_dim <- ", paste0("original_dim <- ", new_visible_dim, " #"), fixed = TRUE, x = tx )
    tx3 <- gsub("intermediate_dim <- ", paste0("intermediate_dim <- ", new_intermediate_dim, " #"), fixed = TRUE, x = tx2 )
    writeLines(tx3, con = file.path(getwd(), "vae_generator_changed.R"))
    warning(paste0("You just changed settings of VAE architecture in '", file,"'. To swap back on old settings you might need to re-source the file."))    
    if (source.file)
        source(file.path(getwd(), "vae_generator_changed.R"))  # it should evaluate at Glob env
}



                                        #

PasteMulticenter <- function(simdlist)
    lapply(1:length(simdlist[[1]]),        ## assuming all list items have same length
           function(h) do.call(
                           "rbind",
                           lapply(1:length(simdlist),
                                  function(tn)
                                      data.frame(simdlist[[tn]][[h]],
                                                 trialID = as.character(tn))
                                  )
                       )
           )



                                        # using all resampling methods
    # simdR is a list of R independent data repetitions

GenerateArtificialDataMultiCenter <- function(simdR, H, site = "trialID",
                                              SI_k = 20000, epochs = 50,
                                              repeat_data = 100, keep = 1:3,
                                              seedAD = 19, DoVae = FALSE, DoVaeOnly = FALSE )
{
    sitenames <- levels(as.factor(simdR[[1]][[site]]))
    R <- length(simdR)
    VAE  <- NA
    if (DoVae)
    {
        VAE <- lapply(1:R,
                      function(r)
                          PasteMulticenter(
                              lapply(sitenames,
                                     function(s)
                                     {   # VAE
                                         D <- simdR[[r]]
                                         set.seed(seedAD + r + which(sitenames == s))
                                         print(
                                             system.time(
                                                 res <- try(
                                                     VAE(D[D[[site]] == s, keep],
                                                         H, epochs, repeat_data ),
                                                     silent = TRUE )
                                             )
                                         )
                                         if (class(res) == "try-error")
                                             return( lapply(1:H,
                                                            function(x)
                                                                replace(simd[1,], values = NA)
                                                            )
                                                    )
                                         else
                                             return(res)
                                     }
                                     )
                          )
                      )
        if (DoVaeOnly)
            return(list(VAE = VAE))
    }

###    IPDboot <- mclapply(1:R, function(r) IPDbootDat(simdR[[r]], H, FALSE, site, seedAD + r) )     # IPD boot (size H), R repetitions

    IPDbootDist <- mclapply(1:R, function(r) IPDbootDat(simdR[[r]], H, TRUE, site, seedAD + r) )      # IPD boot Dist
                                            
    GCgamma <- mclapply(1:R,
                        function(r)
                            PasteMulticenter(
                                mclapply(sitenames,
                                         function(s)
                                         {
                                        # GC gamma
                                             D <- simdR[[r]]
                                             print(
                                                 system.time(
                                                     res <- GCboot(
                                                         D[D[[site]] == s, keep],
                                                         H,
                                                         SI_k,
                                                         3,
                                                         seedAD + r + which(sitenames == s))
                                                 )
                                             )
                                             return(res)
                                         }
                                         )
                            )
                        ) 
                                        
    GCjohn <- mclapply(1:R,
                       function(r)
                           PasteMulticenter(
                               mclapply(sitenames,
                                        function(s)
                                        {   # GC john
                                            D <- simdR[[r]]
                                            print(
                                                system.time(
                                                    res <- GCboot(
                                                        D[D[[site]] == s, keep],
                                                        H,
                                                        SI_k,
                                                        4,
                                                        seedAD + r + which(sitenames == s))
                                                )
                                            )
                                            return(res)
                                        }
                                        )
                           )
                       )
                                        
     
    CART <- mclapply(1:R,
                     function(r)
                         PasteMulticenter(
                             mclapply(sitenames,
                                      function(s)
                                      {   # CART
                                          D <- simdR[[r]]
                                          print(
                                              system.time(
                                                  res <- cart(D[D[[site]] == s, keep],
                                                              H,
                                                              seedAD + r + which(sitenames == s))
                                              )
                                          )
                                          return(res)
                                      }
                                      )
                         )
                     )
    
     IPDboot <- NA### GCgamma  <-  GCjohn  <-  CART <- NA   ### 29.01.2021 deactivating ipdboot for sim_study_v3. MUST reaactivate later .. DELETE
                                        

                                             # vae freezes within mclapply, must use lapply instead
    return( list(IPDboot = IPDboot,
                 IPDbootDist = IPDbootDist,
                 GCgamma = GCgamma,
                 GCjohn = GCjohn,
                 CART = CART,
                 VAE = VAE)
           )
}



## inferences to compute 

Simulations_stat <- function( x, formula,
                             const.haz = FALSE,
                             frailty = FALSE, KM = FALSE,
                             cumhaz = FALSE )
{
    x <- as.data.frame(x)    
    if (!KM){
        if (!frailty){
            if (!const.haz)
                out <- try( coxph( formula , data = x ) , silent = TRUE  )
            else
                out <- try( glm( formula, "poisson", data = x ) , silent = TRUE  )
            if (class(out)[1] == "try-error")
                return( list( coef= NA, sd = NA, vcov = NA, ranefsd= NA, ranefsdSE = NA, revcov = NA, ranef = NA, aic = NA ))
            covmat <- vcov(out)
            return(list(
                coef= coef(out),
                sd = sqrt(diag(covmat)),
                vcov = covmat[lower.tri(covmat)],
                ranefsd= NA,
                ranefsdSE = NA,
                revcov = NA,
                ranef = NA,
                aic = AIC(out) )
                )
        } else {
            out <- try( coxme( formula , data = x ) , silent = TRUE  )
            if (class(out) == "try-error")
                return( list( coef= NA, sd = NA, vcov = NA, ranefsd= NA, ranefsdSE = NA, revcov = NA, ranef = NA, aic = NA ))
            covmat <- vcov(out)
            recovmat <-  as.matrix(VarCorr(out)[[1]])
            return(list(
                coef= fixef(out),
                sd = sqrt(diag(covmat)),
                vcov = covmat[lower.tri(covmat)],
                ranefsd= sqrt(diag(recovmat)),
                ranefsdSE = NA,
                revcov = recovmat[lower.tri(recovmat)],
                ranef = ranef(out)[[1]], aic = AIC(out) )
                )
        }
    }else{
        estim <- ifelse(cumhaz, "cumhaz", "surv")
        out <- try( survfit( formula, x ), silent = TRUE)                 
        if (class(out) == "try-error")
            return(data.frame(km = NA, time = NA, km_l = NA, km_u = NA, trt = c(0,1)))        
        if (!is.null(out$strata)){
            ncontrol <- out$strata[1]
            ntrt <- out$strata[2]
            ntreat <- (ncontrol + 1):(ncontrol + ntrt)
                                        #  to get cumhaz CIs apply -log(surv) on km_l, km_u
            return(rbind(
                data.frame( km = out[[estim]][1:ncontrol], time = out$time[1:ncontrol],
                           km_l = out$lower[1:ncontrol], km_u = out$upper[1:ncontrol], trt = 0 ),
                data.frame( km = out[[estim]][ntreat], time = out$time[ntreat], km_l = out$lower[ntreat],
                           km_u = out$upper[ntreat], trt = 1 )
            )  )
        } else {
            return( data.frame( km = out[[estim]], time = out$time,  km_l = out$lower, km_u = out$upper, trt = NA )  )
        }
    }

}


###

statma <- function(TE, seTE, hubsnames, effect = "HR", statist = c("random", "fixed"))
{
    statist <- match.arg(statist)
    drop <- which(seTE >= 1e+02)
    if (length(drop) > 0)
        out <- try( metagen(TE[-drop], seTE[-drop], hubsnames[-drop], sm = effect, method.tau = "SJ"), silent = TRUE)
    else
        out <- try( metagen(TE, seTE, hubsnames, sm = effect, method.tau = "SJ"), silent = TRUE)
    if (class(out)[[1]] == "try-error")
       return( data.frame( coef= NA, sd = NA, vcov = NA, ranefsd= NA, ranefsdSE = NA, revcov = NA, ranef = NA, aic = NA ))
    else
        return( data.frame(
            coef= out[[paste0("TE.",statist)]],
            sd = out[[paste0("seTE.",statist)]],
            ranefsd = ifelse(statist[1] == "random", out$tau, NA),
            ranefsdSE = ifelse(statist[1] == "random", out$se.tau2, NA),
            ranefcor = NA, aic = NA ) )       
}





CollectMAsummaries <- function(simd, formula = NULL, site = "trialID", noMA = FALSE)
{
    out <- NULL
    if (!is.null(formula))
    {
        sitenames <- levels(as.factor(simd[[site]]))
        out <- lapply(sitenames, function(s){
            Simulations_stat(simd[simd[[site]] == s, -which(names(simd) == site)], formula)
        }
        )
        names(out) <- sitenames
        MAdat <- lapply(names(out[[1]]$coef),
                        function(b) do.call("rbind",
                                            lapply(sitenames, function(s)
                                                data.frame( beta = out[[s]]$coef[b],
                                                           betaSD = out[[s]]$sd[b],
                                                           study = s)
                                                )
                                            )
                        )
        names(MAdat) <- names(out[[1]]$coef)
    }
    if (is.null(out)){
        sitenames <- names(simd)
        out <- simd
    }
    if (noMA)
        MAdat <- do.call(
            "rbind",
            lapply(c("coef", "sd"), function(i)
                do.call(
                    "rbind",
                    lapply(names(out[[1]]$coef), function(b)
                        do.call(
                            "rbind",
                            lapply(sitenames, function(s)
                                data.frame(est = switch(i,
                                                         "coef" = "$\\hat\\beta$",
                                                        "sd" = "$\\hat\\sigma_{\\hat\\beta}$"
                                                        ),
                                           variab = b,
                                           site = s,
                                           ipd.ref =  out[[s]][[i]][b]
                                           )
                                )
                        )
                        )
                )
                )
        )
    return(MAdat)
}
    
                                        #

MultiplMA <- function(MAdat, MAeffect, statist)
{
    MAinfer <- lapply(MAdat, function(x) statma(x$beta, x$betaSD, x$study, MAeffect, statist))
    variab <- names(MAdat)
    params <- names(MAinfer[[1]])
    res <- lapply(params, function(b)
    {
        out <- unlist(lapply(MAinfer, function(x) x[[b]]) )
        names(out) <- variab
        return(out)
    }
    )
    names(res) <- params
    return(res)
}

                                        #

MA <- function(simd, formula = NULL, site = "trialID",
               MAeffect = "HR", statist = c("random", "fixed") )
    MultiplMA(CollectMAsummaries(simd, formula, site), MAeffect, statist)

                                        #


 ##                                        #  it works both on IPD or on pseudodata

Simulations_BootStat <- function(bootdat, statist = c("random", "consthaz", "km"),
                                 reg_formula = NULL, trim.ol = FALSE, cumhaz = FALSE )
{
    statist <- match.arg(statist)
                                        # bootstrap statistics
    IPDboot <-  mclapply(bootdat, function(x)
        Simulations_stat(x, reg_formula,
                         ifelse(statist == "consthaz", T, F ),
                         ifelse(statist == "random", T, F ),
                         ifelse(statist == "km", T, F),
                         cumhaz)
        )
    if (statist != "km")
    {
        Bag <- BagInfer(IPDboot, "mean", trim.ol)
        BagSE <- BagInfer(IPDboot, "sd", trim.ol)
        Bag95QQ <- BagInfer(IPDboot, "95qq", trim.ol)
    }
    else
    {
        Bag <- AverageKM(IPDboot)
        BagSE <- NA ## here probably prone to errors later
        Bag95QQ <- NA ## here probably prone to errors later
    }
    return( list( Bag = Bag, BagSE = BagSE, Bag95QQ = Bag95QQ) )
}


## SIMULATION MODULE (data + inference)   ## changed 24.6.20 arg 'rdn_sd'
     
Simulations <- function(R = 1, K = 2, n = 100, gam, lam, beta, rdn_sd = c(0.7, 0.0),    ### rho upcoming feat
                        rho = NULL, rateC = c(0.03, 0.05), second_scenario = FALSE, 
                        reg_formula = Surv( time, status )~( 1|as.factor(trialID) ) + trt,
                        statist = c("random", "consthaz", "km"), H= 100, SI_k = 20000, epochs = 50,
                        repeat_data = 100, trim.ol = FALSE, keep = 1:3, seed = 51,
                        DoVae = FALSE, DoVaeOnly = FALSE, noRefan = FALSE)
{
    statist <- match.arg(statist)
    fixedform <- reg_formula
    fixedform[3] <- expression(trt)        
    ## simulate IPD (R independent trials)
    simdR <- generateDataScenario(R, K, n, gam, lam, beta, rdn_sd, rho, rateC, second_scenario, seed)                                                               # generate artificial data (HEAVY COMPUTATIONS HERE)
    ad <-  GenerateArtificialDataMultiCenter(simdR,
                                             H,
                                             "trialID",
                                             SI_k,
                                             epochs,
                                             repeat_data,
                                             keep,
                                             seed,
                                             DoVae,
                                             DoVaeOnly)

                                        # compute inferences and average
    
    ### saveRDS(ad, "simul_3_VAE_raw_.rds")   ### interim saving DELETE
    
    Coxipd <- Coxma <- list(MC_MEAN = NA, MC_SD = NA)
    KMipd  <- KMma <- data.frame(km = NA, time = NA, km_l = NA, km_u = NA, trt = NA)
    if (!noRefan)
    { 
        Coxipds <- mclapply( simdR, function(d)
            Simulations_stat(d, reg_formula,
                             ifelse(statist == "consthaz", T, F),
                             ifelse(statist == "random", T, F),
                             FALSE
                             )
            )
        Coxipd <- list(MC_MEAN = BagInfer(Coxipds), MC_SD = BagInfer(Coxipds, "sd"))
        KMipd <- AverageKM(mclapply( simdR, function(d) Simulations_stat(d, fixedform, FALSE, FALSE, TRUE ) ) )
        Coxmas <- lapply(simdR, function(d) MA(d, fixedform)) 
                                        #      loghrs <- lapply(Coxmas, function(x) x$coef)  # use loghrs in param KM below (DEACTIVATED)
        Coxma <- list( MC_MEAN = BagInfer(Coxmas), MC_SD = BagInfer(Coxmas, "sd") )
        KMma <- AverageKM( lapply(1:R, function(d)
            kmparam(simdR[[d]][["time"]],
                    simdR[[d]][["status"]],
                    simdR[[d]][["trt"]], NULL,
                    simdR[[d]][["trialID"]]
                    )
            )
            )  # parametric KM (no random, no loghr)
    }
    Coxreg <- lapply(ad, function(m)
    {
        res <- mclapply(m, function(d)
            Simulations_BootStat(d, statist, reg_formula)
            )  # here cannot be "km"
                                        # return experimental averages and SDs
        return( list(
            MC_MEAN = list(Bag = BagInfer(lapply(res, function(x) x$Bag)),
                           BagSE = BagInfer(lapply(res, function(x) x$BagSE)),
                           Bag95QQ = BagInfer(lapply(res, function(x) x$Bag95QQ)) ),
            MC_SD = list(Bag = BagInfer(lapply(res, function(x) x$Bag), "sd"),
                         BagSE = BagInfer(lapply(res, function(x) x$BagSE), "sd"),
                         Bag95QQ = BagInfer(lapply(res, function(x) x$Bag95QQ), "sd" )
                         )
        )
        )
    }
    )
    KMcurves <- lapply(ad, function(m) AverageKM(mclapply(m, function(d)
        Simulations_BootStat(d, "km", fixedform)$Bag
        )
        )
        )
    rm(ad, simdR)  # memory clean-up !!!
    return( list(Coxipd = Coxipd, KMipd = KMipd,
                 Coxreg = Coxreg, KMcurves = KMcurves,
                 Coxma = Coxma, KMma = KMma) )
}


                                        #
TableReg <- function(regres, colname)
{
    names <- c("coef", "sd", "ranefsd")
    out <- data.frame("MC_MEAN" = unlist(regres$MC_MEAN[names]),
                      "MC_SD" = unlist(regres$MC_SD[names])
                      )
    colnames(out) <- paste(colname, colnames(out), sep = ".")
    return(out)
}

 
                                        #

TableRes <- function(simres, moreinfo = NULL)
{
    KM <- rbind(data.frame(simres$KMipd,
                           type = "Orig. IPD"),
                do.call(
                    "rbind",
                    lapply(names(simres$KMcurves), function(m)
                        data.frame(simres$KMcurves[[m]],
                                   type = m
                                   )
                        )
                ),
                data.frame(simres$KMma,
                           type = "Meta an."
                           )
                )
    REG <- cbind(data.frame(
        Est = c("coef", "sd", "ranefsd"),
        TableReg(simres$Coxipd, "IPD")
    ),
    do.call(
        "cbind",
        lapply(names(simres$Coxreg), function(m)
            TableReg(lapply(simres$Coxreg[[m]],
                            function(x) x$Bag), m)
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



 ### TODO 24.6.20 check new version above for compatibility with old one below 
## generateDataScenario <- function(R = 1, K = 2, n = 1000, gam = c(1.5, 1.05), lam = c(0.1, 0.27),
##                                  beta = c(-0.5, -0.005), rdn_sd = 0.7, rateC = c(0.03, 0.05),
##                                  second_scenario = FALSE, seed = 61 ){

##     n <- ifelse(second_scenario, n/2, n)
    
##   simd <- mclapply(1:R, function(i){   #### list of R independent experiments 

##    set.seed(seed + i)
      
## frailty <- rgmom(K, mean = 1, sd = rdn_sd) # frailty distr

##       s <- do.call("rbind", mclapply(frailty, function(j){
          
##  set.seed(seed + i + j)
##  s1 <- data.frame( Weibtimes(n, lam[1], gam[1], beta[1], rateC[1], j ), trialID = which(frailty == j), scen = 1 )
##  s2 <- NULL
##  if (second_scenario)
##  s2 <- data.frame( Weibtimes(n, lam[2], gam[2], beta[2], rateC[2], j ), trialID = which(frailty == j), scen = 2 )
      
##   return(rbind(s1, s2))
 
##       }   ) )
##     }   )

##           return(simd)  # list of R simulated datasets       
##   }





## vae <- function(simd, H, epochs = 50, repeat_data = 100, seedAD){

##    source("vae_generator_changed.R", local = TRUE)   # BEWARE: u must source in-loop. Apparently local envir forgets GlobEn compiling (?)
##     set.seed(seedAD)
    
##     out  <- try( VAE(simd, H, epochs, repeat_data), silent = TRUE )

##               if (class(out) == "try-error")
##         return( lapply(1:H, function(x) replace(simd[1,], values = NA) )  ) 
    
##     return(out)
##     }










