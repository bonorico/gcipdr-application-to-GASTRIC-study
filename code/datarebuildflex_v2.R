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




#' @title Tailored GC optimization
#' @description Sequential GC optimization, based on analytic Kruskal conversion, or alternating integration (stochastic or numerical).
#' @param SI_k_seq a sequence of (increasing) values for argument 'SI_k' of DataRebuild.
#' @param NI_maxEval_seq a sequence of (increasing) values for argument 'NI_maxEval' of DataRebuild.
#' @param plan_switch integer atomic vector - If equals 1, it indicates numerical integration. If 2, it indicates stochastic integration.
#'    E.g. c(1,2,2,1) defines the plan = numeric, stochastic, stochastic, numeric, in this specific order. The length of 'NI_maxEval_seq' ('SI_k_seq
#'    ') must agree with number of 1's (2's).
#' @param try_kruskal_first logical - Should Kruskal conversion be performed ? Default FALSE. If TRUE it overrides all other integration options.
#' @param treshold numerical - Threshold in absolute value for the accepted difference between observed and simulated correlation. Default = 0.15. 
#' @details For all other arguments refer to 'DataRebuild'.
#' @seealso DataRebuild

                                        # Test 26.12.2021

tailored_optim_GC <- function(H, n, correlation.matrix, moments, x.mode,
                              corrtype = c("moment.corr", "rank.corr", "normal.corr"),
                              marg.model = c("gamma", "johnson"), variable.names = NULL,
                              SBjohn.correction = FALSE, compute.eec = FALSE, checkdata = FALSE,
                              tabulate.similar.data = FALSE, SI_k_seq = c(8000, 30000, 50000),
                              NI_tol = 0.01, NI_maxEval_seq = c(500, 2000), plan_switch = c(1, 2, 2, 1, 2),
                              try_kruskal_first = FALSE, treshold = 0.15, s.seed = NULL)
{
    if (length(plan_switch[plan_switch == 1]) != length(NI_maxEval_seq) |
        length(plan_switch[plan_switch == 2]) != length(SI_k_seq))
        stop("in plan_switch: Number of '1's ('2's) must equal length of 'NI_maxEval_seq' ('SI_k_seq').")
    performance <- c()  
    logmex <- c()  # log message
    stop_condition <- quote(ifelse(class(out) == "similar.data",
                                       all(abs(out$is.similar[[5]][, 3]) <= treshold), FALSE) )  # threshold between in and out correlation
    record_performance <- quote( ifelse(class(out) == "similar.data", 
                                        mean(abs(out$is.similar[[5]][, 3]) > treshold), NA ) )  # record mean failure incidence  
    mex <- paste0( c("Kruskal", "Kruskal fine-tuned", "Numerical", "Stochastic"), " integration. ")
    names(mex) <- as.character(-1:2)
    plan_switch <- c(-1:0, plan_switch)
    plan_N <- length(plan_switch)
    optim_seq <- vector(mode = "integer", length = plan_N)
    optim_seq[which(plan_switch == 1)] <-  NI_maxEval_seq
    optim_seq[which(plan_switch == 2)] <-  SI_k_seq
    for (i in 1:plan_N)
    {
        if (!try_kruskal_first & plan_switch[i] < 1)     # skip if Kruskal approach is set to FALSE
            next
        planindex <- as.character(plan_switch[i])
        current_plan <- paste0(mex[planindex],
                               c("1" = "NI_maxeval = ", "2" = "SI_k = ")[planindex], optim_seq[i] )
        print(current_plan)
        logmex <- c(logmex, current_plan)
        set.seed( s.seed, "L'Ecuyer")
        out  <- try(
            DataRebuild(H, n, correlation.matrix, moments,
                        x.mode, marg.model = marg.model,
                        variable.names = variable.names, checkdata = checkdata,
                        tabulate.similar.data = tabulate.similar.data,
                        stochastic.integration = ifelse(plan_switch[i] == 2, TRUE, FALSE),
                        SI_k = optim_seq[i], NI_maxEval = optim_seq[i],
                        assume.all.smooth = ifelse(try_kruskal_first & plan_switch[i] == -1, TRUE, FALSE),
                        cp.finetune = ifelse(try_kruskal_first & plan_switch[i] == 0, TRUE, FALSE)),
            silent = TRUE )
        performance <- c(performance, eval(record_performance))
        if ( eval(stop_condition) )
            break
    }
    print(paste0("Exiting with ", logmex[length(logmex)]))
    return(list(res = out,
                logmex = data.frame(exitmex = logmex[length(logmex)],
                                    perfexit = performance[length(performance)],
                                    bestperf = logmex[which(performance == min(performance, na.rm = TRUE))[1]],
                                    bestperfval = min(performance, na.rm = T)
                                    )
                )
           )

}  




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


