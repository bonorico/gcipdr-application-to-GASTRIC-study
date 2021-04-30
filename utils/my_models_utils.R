## my_models_utils.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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




####

check <- function(x)
{
    d <- which(is.null(x) | length(x) < 1)
    x[d] <- NA
    return(x)
}



#### changed 01.12.2020 add Schoenfeld residuals stat

Schoenfeld_resid_and_test <- function(dat, formula)
{

    out <- try(
        cox.zph(
            coxph(formula, dat, x = TRUE),
            global = F
        )
    )
    if (class(out)[1] == "try-error")
        return(list(chisq = NA, pval = NA))  ## compute p-val directly on averaged chisq
    else
        return(list(chisq = out$table[, 1], pval = NA))

  }

                                        #  IMPORTANT PENDING --> edit code in _v13 'Simulations' OR separate objects ..

stat <- function( x, formula,
                 model = c("coxme", "coxph", "cox.zph", "glm"),
                 KM = FALSE, cumhaz = FALSE )
{

    model <- match.arg(model)

    x <- as.data.frame(x)
    
    if (!KM)
    {
        if (model == "cox.zph")
            return(Schoenfeld_resid_and_test(x, formula))  ### just Schoenfeld test on coxph residuals
        
        out <- switch( model,
                      "glm" = try( glm( formula, "poisson", data = x ) , silent = TRUE  ),
                      "coxph" = try( coxph( formula , data = x ) , silent = TRUE  ),
                      "coxme" = try( coxme( formula , data = x ) , silent = TRUE  )  )
        if (class(out)[1] == "try-error")
            return( list(coef= NA,
                         sd = NA,
                         vcov = NA,
                         ranefsd2= NA,
                         ranefsd1= NA,
                         ranefsdSE = NA,
                         revcov = NA,
                         ranef = NA,
                         lr = NA,
                         df = NA
                         )
                   )
        covmat <- vcov(out)       

        revcov <-  NA
        ranefsd1 <-  NA
        ranefsd2 <-  NA
        ranef <- NA
        
        if (model == "glm")
        {
            lr <- out$deviance
            df <- out$df.residual
        }
        if (model == "coxph")
        {
            sout <- summary(out)
            lr <- sout$logtest[["test"]]
            df <- sout$logtest[["df"]]
        }   
        if (model == "coxme")
        {
            lr <- 2*(out$loglik[["Integrated"]] - out$loglik[["NULL"]])
            df <- out$df[1]
            recovmat <- as.matrix(VarCorr(out)[[1]])
            revcov <- recovmat[lower.tri(recovmat)]
            ranefsd <-  sqrt(diag(recovmat))
            ranefsd1 <-  ranefsd[1]
            ranefsd2 <-  ifelse(length(ranefsd) > 1, ranefsd[2], NA)
            ranef  <-  ranef(out)[[1]]
        }


        res <- list(coef= coef(out),
                    sd = sqrt(diag(covmat)),
                    vcov = covmat[lower.tri(covmat)],
                    ranefsd2= ranefsd2,
                    ranefsd1= ranefsd1,
                    ranefsdSE = NA,
                    revcov = revcov,
                    ranef = ranef,
                    lr = lr,
                    df = df
                    )
        
        return(
            lapply(res, function(i) check(i))
        )
    }
    else
    {
        estim <- ifelse(cumhaz, "cumhaz", "surv")
        out <- try(
            survfit( formula, x ),
            silent = TRUE
        )                 

        if (class(out) == "try-error")
            return(data.frame(km = NA,
                              time = NA,
                              km_l = NA,
                              km_u = NA,
                              trt = c(0,1)
                              )
                   )
        
        if (!is.null(out$strata))
        {
            ncontrol <- out$strata[1]
            ntrt <- out$strata[2]
            ntreat <- (ncontrol + 1):(ncontrol + ntrt)
                                        #  to get cumhaz CIs apply -log(surv) on km_l, km_u
            return(
                rbind(
                    data.frame(km = out[[estim]][1:ncontrol],
                               time = out$time[1:ncontrol],
                               km_l = out$lower[1:ncontrol],
                               km_u = out$upper[1:ncontrol], trt = 0 ),
                    data.frame(km = out[[estim]][ntreat],
                               time = out$time[ntreat],
                               km_l = out$lower[ntreat],
                               km_u = out$upper[ntreat], trt = 1
                               )
                )
            )
        }
        else
        {
            return(
                data.frame(km = out[[estim]],
                           time = out$time,
                           km_l = out$lower,
                           km_u = out$upper,
                           trt = NA
                           )
            )
        }

    }

}


##


## CHANGED 01.12.2020  adding cox.zph                                        #  it works both on IPD or on pseudodata

BootStat <- function(bootdat,
                     statist = c("random", "fixed", "residuals", "consthaz", "km"),
                     reg_formula = NULL, trim.ol = FALSE, cumhaz = FALSE )
{
    statist <- match.arg(statist)
                                        # bootstrap statistics
    IPDboot <-  mclapply(bootdat, function(x)
        stat(x, reg_formula, ifelse(statist == "random", "coxme",
                             ifelse(statist == "fixed", "coxph",
                             ifelse(statist == "residuals", "cox.zph", "glm") ) ),  # overrun by arg "km"
             ifelse(statist == "km", T, F), cumhaz)
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
    return( list(Bag = Bag,
                 BagSE = BagSE,
                 Bag95QQ = Bag95QQ
                 )
           )

}




###  AVERAGE KM (NA) estimates

AverageKM <- function(MClist)
{
    out <-  try(
        Average.MC.cumhaz(MClist),
        silent = TRUE
    )
    if (class(out) == "try-error")
        return(data.frame(km = NA,
                          time = NA,
                          km_l = NA,
                          km_u = NA,
                          trt = c(0,1)
                          )
               )
    else
        return(out$mc.average)
}



## PARAMETRIC surv func extimate

survpp <- function(t, alpha, alphaSD)
{
    St <- exp(-alpha*t)
    
    z <- qnorm(0.025, lower.tail = FALSE)
    
    logSD <- sqrt( ((-t)^2)*alphaSD^2 )  # delta method on log exp(-alpha*t) wrt alpha

    logLow <- log(St) - (z*logSD)
    logUp <- log(St) + (z*logSD) 

    return( list(St = St,
                 low = exp(logLow),
                 up = exp(logUp)
                 )
           )
}  
                                        


mle_haz <- function(status, time)
{
    rs <- sum(time)

    alpha <- sum(status)/rs

    alpha_sd <- alpha/rs

    return( data.frame(haz = alpha,
                       hazSD = alpha_sd
                       )
           )
    
}
                                        


MLEhaz <- function(status, time, site = NULL)
{
    if (!is.null(site))
    {
        dat <- do.call(
            "rbind",
            lapply(levels(as.factor(site)), function(s)
                data.frame(events = sum(status[site == s]),
                           rs = sum(time[site == s]),
                           site = s)
                )
        )
        out <-  metarate(events, rs, site, dat, sm = "IR", method.tau = "SJ")
        fixhaz <- out$TE.fixed
        fixhazSD <- out$seTE.fixed
        mixhaz <- out$TE.random
        mixhazSD <- out$seTE.random
    }
    else
    {
        haz <- mle_haz(status, time)
        fixhaz <- haz$haz
        fixhazSD <- haz$hazSD
        mixhaz <- NA
        mixhazSD <- NA
    }

    return( list(fixhaz = fixhaz,
                 fixhazSD = fixhazSD,
                 mixhaz = mixhaz,
                 mixhazSD = mixhazSD
                 )
           )
}


  #
kmparam <- function(time, status,
                    trt, maxt = NULL,
                    site = NULL,
                    random.effect = FALSE,
                    loghr = NULL,
                    loghrSD = NULL)
{
    study0 <- study1 <- NULL
    if (!is.null(site))
    {
        study0 <- site[trt == 0]
        study1 <- site[trt == 1]
    }
    
    out0 <- MLEhaz(status[trt == 0], time[trt == 0], study0)
    out1 <- MLEhaz(status[trt == 1], time[trt == 1], study1)

    if (random.effect)
    {
        haz0 <- out0$mixhaz
        haz0SD <- out0$mixhazSD 
        haz1 <- out1$mixhaz
        haz1SD <- out1$mixhazSD
    }
    else
    {
        haz0 <- out0$fixhaz
        haz0SD <- out0$fixhazSD
        haz1 <- out1$fixhaz
        haz1SD <- out1$fixhazSD
    }

    if (!is.null(loghr) & !is.null(loghrSD))
    {
        haz1 <- haz0*exp(loghr)
        haz1SD <- exp(loghrSD) # delta method on log(haz1) wrt loghr
     }

    tt <- seq(0.01,
              ifelse(is.null(maxt), max(time, na.rm = T), maxt),
              length = 1000)

    St0 <- survpp(t = tt, alpha = haz0, alphaSD = haz0SD )
    St1 <- survpp(t = tt , alpha = haz1, alphaSD = haz1SD)

    return(rbind(
        data.frame(
            km = St0$St ,
            time = tt,
            km_l = St0$low,
            km_u = St0$up,
            trt = 0 ),
        data.frame(km = St1$St,
                   time = tt,
                   km_l = St1$low,
                   km_u = St1$up, trt = 1
                   )
    )
    )  

}
##



### BOOT UTILS

#### BAGGING ESTIMATES

BagInfer <- function(inferencesSample,
                     operator = c("mean", "sd", "95qq"),
                     trim.ol = FALSE )
{
    operator <- match.arg(operator)
    op <- function(type, x) switch(type,
                                   "mean" = mean(x, na.rm = TRUE),
                                   "sd" = sd(x, na.rm = TRUE),
                                   "95qq" = quantile(x, c(0.025, 0.975), na.rm = TRUE) )
     
    betas <- names(inferencesSample[[1]])
    out <- lapply(betas, function(b)
        apply(do.call(
            "rbind",
            lapply(inferencesSample, function(y) y[[b]]) ), 2,
            function(x) op(operator, trimOutliers(x, trim.ol))))
    names(out) <- betas

    return(out)
    
}


#### eventually trim simulation outliers (typically set to FALSE)

trimOutliers <- function(x, trim = TRUE)
{
                                        # discard outliers
    if (!trim)
        return(x)
    m <- median(x,na.rm=T)
    s <- mad(x, na.rm=T)
        
    if (all(is.na(m)) | all(is.na(s)))
    {
        warning("discard: cannot compute median or mad (all x values missing)")
        return(x)
     }
        
    if(max(x, na.rm=T) > (m+s*3.5) | min(x,na.rm=T) < (m-s*3.5))
    {    
    min <- m-(2.5*s)
    max <- m+(2.5*s)

    res <- x[x>=min & x<=max]

    }else
    {
        res <- x
    }
    res
}



## surrogates endpoints 

surro <- function(d)
{
    inp <- data.frame(id = 1:dim(d)[1],
                      trialref = d[, "trialID"],
                      d[, -which(names(d) == "trialID")]
                      )
    inp$timeT <- ifelse(inp$timeT < 0, 0.0001, inp$timeT )
    inp$timeS <- ifelse(inp$timeS < 0, 0.0001, inp$timeS )

    out <- try(
        surrosurv(inp, model = "clayton"),
        silent = TRUE)

    if (class(out)[1] == "try-error")
        return(list(theta = NA, R2 = NA, logSTE = NA, conv = NA))
    
    theta <- out[[1]]$unadj$step1$theta  ## individual surrogacy
    R2 <- out[[1]]$unadj$R2   # trial-level surrogacy
    logSTE <- surrosurv::ste(out)[["Clayton.unadj"]]  # surrogate treshold effect (lower pred interv for HR true endpoint)
    conv = all(convergence(out)[1,])

    return(list(theta = theta, R2 = R2, logSTE = logSTE, conv = conv))    

}

#


surroboot <- function(d, trim.ol = FALSE)
{                                        # bootstrap statistics
    out <- mclapply(d, function(x)
    {
        res <- try(
            surro(x),
            silent = F)  ### must re-try to by-pass internal settings ...
        if (class(res)[1] == "try-error")
            return(list(theta = NA, R2 = NA, logSTE = NA, conv = NA))
        else
            return(res)      
    }
    )    

    Bag <- BagInfer(out, "mean", trim.ol)
    BagSE <- BagInfer(out, "sd", trim.ol)
    Bag95QQ <- BagInfer(out, "95qq", trim.ol)
    return( list(Bag = Bag,
                 BagSE = BagSE,
                 Bag95QQ = Bag95QQ)
           )

}






## (IPD) RESAMPLING METHODS
## 27/12/19 changed: deleting argument 'avoid.zero.sd'

IPDbootDat <- function(simd, H, distributed = FALSE, site = "trialID", seedAD = 82 )
{
    set.seed(seedAD)
    if (distributed){  # distributed boot data
        sitenames <- levels(as.factor(simd[[site]]))
        out <- lapply(sitenames, function(s) lapply(1:H, function(h){
            res <-  simd[simd[[site]] == s, -which(names(simd) == site)]
            resboot <- res[sample(1:dim(res)[1], dim(res)[1], TRUE), ]
            return(resboot)
        }
        )
        )
        bootdat <- lapply(1:H, function(h) do.call("rbind", # merge by site
                                                   lapply(1:length(sitenames), function(s){
                                                       bd <- as.data.frame(out[[s]][[h]])  # matrix simul
                                                       bd <- data.frame(bd, sitenames[s])
                                                       colnames(bd)[dim(bd)[2]] <- site
                                                       return(bd)
                                                   }
                                                   )
                                                   )
                          )
    }
    else  # pooled bood data
        bootdat <- lapply(1:H, function(h) simd[sample(1:dim(simd)[1], dim(simd)[1], TRUE), ]   )
    return(bootdat)
    
}
