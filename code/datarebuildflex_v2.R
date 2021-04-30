## datarebuildflex_v2.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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




#' @name PasteDistributed
#' @title Pool center-specific artificial data
#'
#' @description Performs a rudimental adaptive NORTA optimizationg based on rough Kruskal conversion, or stochastic or numerical integration.
#'
#' @param SI_k_seq a sequence of three (increasing) values for the 'SI_k' argument of DataRebuild.
#' 

PasteDistributed <- function(obj, hubsvariable, hubsnames)
{
    if( class(obj[[1]]) == "similar.data" )
        H <- length(obj[[1]][["Xspace"]])
    listofpooled <-  lapply(1:H, function(h)
    {      # by simulation
        merged <- do.call(
            "rbind",
            lapply(1:length(hubsnames), function(s)
            {  ##  bind by country
                ad <- as.data.frame(obj[[s]][["Xspace"]][[h]])  # matrix simul
                ad <- data.frame(ad, hubsnames[s])
                colnames(ad)[dim(ad)[2]] <- hubsvariable
                return(ad)
            }
            )
        )
        return(merged)
    }
    )
    return(listofpooled)
}



#' @name DataRebuildFlex
#' @title Flexible DataRebuild
#'
#' @description Performs a rudimental adaptive NORTA optimizationg based on rough Kruskal conversion, or stochastic or numerical integration.
#'
#' @param SI_k_seq a sequence of three (increasing) values for the 'SI_k' argument of DataRebuild.
#'
#' @param NI_maxEval_seq a sequence of two (increasing) values for the 'NI_maxEval' argument of DataRebuild.
#'
#' @param assume.all.smooth logical. Should rough Kruskal conversion (with fine-tuning) be performed ? Default TRUE.


## adaptive data rebuild based on simulation outcome

DataRebuildFlex <- function(H, n, correlation.matrix, moments, x.mode,
                            corrtype = c("moment.corr", "rank.corr", "normal.corr"),
                            marg.model = c("gamma", "johnson"), variable.names = NULL,
                            SBjohn.correction = FALSE, compute.eec = FALSE, checkdata = FALSE,
                            tabulate.similar.data = FALSE, SI_k_seq = c(8000, 30000, 50000),
                            NI_tol = 0.01, NI_maxEval_seq = c(500, 2000),
                            assume.all.smooth = FALSE, treshold = 0.15, s.seed = NULL)
{
    perf <- c()  # performance
    logmex <- c()  # log message
    out <- NULL
    condition_refused <- quote(
        class(out) == "try-error" |
        ifelse(class(out) == "similar.data",
               any(abs(out$is.similar[[5]][, 3]) > treshold),
               FALSE)
    )    
### mixed approach
    if (assume.all.smooth)
    {
        mex1 <- paste("using Kruskal analytic solution")
        logmex <- c(logmex, paste0("exiting after ", mex1))
        print(mex1)
        out  <- try(
            DataRebuild(H, n, correlation.matrix, moments,
                        x.mode, marg.model = marg.model,
                        corrtype = corrtype, variable.names = variable.names,
                        checkdata = checkdata, tabulate.similar.data = tabulate.similar.data,
                        assume.all.smooth = TRUE ),
            silent = TRUE )
        perf <- c(perf,
                  ifelse(class(out) != "try-error",
                         mean(abs(out$is.similar[[5]][, 3]) > treshold),
                         NA
                         )
                  )
        if ( eval(condition_refused) )
        {
            mex2 <- paste("fine-tuning Kruskal analytic solution")
            logmex <- c(logmex, paste0("exiting after ",mex2))
            print(mex2)
            set.seed( s.seed, "L'Ecuyer")
            out  <- try(
                DataRebuild(H, n, correlation.matrix, moments,
                            x.mode, marg.model = marg.model,
                            variable.names = variable.names, checkdata = checkdata,
                            tabulate.similar.data = tabulate.similar.data,
                            assume.all.smooth = TRUE, cp.finetune = TRUE  ),
                silent = TRUE )
            perf <- c(perf, ifelse(class(out) != "try-error", mean(abs(out$is.similar[[5]][, 3]) > treshold), NA ))

   
                  }
    }
    if (ifelse(!is.null(out),
               eval(condition_refused),
               TRUE)
        )
    {
        mex3 <- ifelse( !is.null(out),
                       paste("switching to numerical integration. NI_maxEval =", NI_maxEval_seq[1]),
                       paste("using numerical integration. NI_maxEval =", NI_maxEval_seq[1]) )
        logmex <- c(logmex, paste0("exiting after ",mex3))
        print(mex3)
        set.seed( s.seed, "L'Ecuyer")
        out  <- try(
            DataRebuild(H, n, correlation.matrix, moments, x.mode,
                        marg.model = marg.model, NI_maxEval = NI_maxEval_seq[1],
                        variable.names = variable.names, checkdata = checkdata,
                        tabulate.similar.data = tabulate.similar.data  ),
            silent = TRUE )
        perf <- c(perf, ifelse(class(out) != "try-error", mean(abs(out$is.similar[[5]][, 3]) > treshold), NA ))
    }
    if ( eval(condition_refused) )
    {    # increase maxevals
        mex4 <- paste("switching to stochastic integration. SI_k =", SI_k_seq[1])
        logmex <- c(logmex, paste0("exiting after ",mex4))
        print(mex4)
        set.seed( s.seed, "L'Ecuyer")
        out  <- try(
            DataRebuild(H, n, correlation.matrix, moments,
                        x.mode, marg.model = marg.model, stochastic.integration = TRUE,
                        SI_k = SI_k_seq[1], variable.names = variable.names,
                        checkdata = checkdata, tabulate.similar.data = tabulate.similar.data  ),
            silent = TRUE )
        perf <- c(perf, ifelse(class(out) != "try-error", mean(abs(out$is.similar[[5]][, 3]) > treshold), NA ))
    }
    if ( eval(condition_refused) )
    {    # increase maxevals
        mex5 <- paste("increasing SI_k =", SI_k_seq[2])
        logmex <- c(logmex, paste0("exiting after ",mex5))
        print(mex5)
        set.seed( s.seed, "L'Ecuyer")
        out  <- try(
            DataRebuild(H, n, correlation.matrix, moments, x.mode,
                        marg.model = marg.model, stochastic.integration = TRUE,
                        SI_k = SI_k_seq[2], variable.names = variable.names,
                        checkdata = checkdata, tabulate.similar.data = tabulate.similar.data  ),
            silent = TRUE )
        perf <- c(perf, ifelse(class(out) != "try-error", mean(abs(out$is.similar[[5]][, 3]) > treshold), NA ))

    }
    if ( eval(condition_refused) )
    {   # increase maxevals
        mex6 <- paste("switching to numerical integration. NI_maxEval =", NI_maxEval_seq[2])
        logmex <- c(logmex, paste0("exiting after ",mex6))
        print(mex6)
        set.seed( s.seed, "L'Ecuyer")
        out  <- try(
            DataRebuild(H, n, correlation.matrix, moments, x.mode,
                        marg.model = marg.model, NI_maxEval = NI_maxEval_seq[2],
                        variable.names = variable.names, checkdata = checkdata,
                        tabulate.similar.data = tabulate.similar.data  ),
            silent = TRUE )
        perf <- c(perf, ifelse(class(out) != "try-error", mean(abs(out$is.similar[[5]][, 3]) > treshold), NA ))

    }
    if ( eval(condition_refused) )
    {    # increase maxevals
        mex7 <- paste("switching to stochastic integration. SI_k =", SI_k_seq[3])
        logmex <- c(logmex, paste0("exiting after ",mex7))
        print(mex7)
        set.seed( s.seed, "L'Ecuyer")
        out  <- try(
            DataRebuild(H, n, correlation.matrix, moments, x.mode,
                        marg.model = marg.model, stochastic.integration = TRUE,
                        SI_k = SI_k_seq[3], variable.names = variable.names,
                        checkdata = checkdata, tabulate.similar.data = tabulate.similar.data),
            silent = TRUE )
        perf <- c(perf, ifelse(class(out) != "try-error", mean(abs(out$is.similar[[5]][, 3]) > treshold), NA ))
        
    }
    print(logmex[length(logmex)])
    return(list(res = out,
                logmex = data.frame(exitmex = logmex[length(logmex)],
                                    perfexit = perf[length(perf)],
                                    bestperf = logmex[which(perf == min(perf, na.rm = T))[1]],
                                    bestperfval = min(perf, na.rm = T)
                                    )
                )
           ) 
    
}





