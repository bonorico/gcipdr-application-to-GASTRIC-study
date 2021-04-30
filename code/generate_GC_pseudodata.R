## generate_GC_pseudodata.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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




######  PSEUDODATA generation (copula simulation)



## PHASE 2: generate artificial IPD from IPD summaries, using copula


GC_pseudodata <- function()
{
    H <- 300  # set simulation size
    seed <- 734   #  734932
    system.time(
                                        # NOTE: can just run once --> next time load saved artificial data (below) directly
        artds <-  lapply(
            c("gamma", "johnson"), function(j)
            {
                out <- lapply(trialnames, function(tn)
                {
                    s <-  trial_summaries[[tn]] # trial data
                    jiseed <- as.integer(tn) + seed + nchar(j)
                    set.seed( jiseed, "L'Ecuyer") # delete this line to assess stability
                    print(
                        system.time(
                            artificial.data <- DataRebuildFlex(H,
                                                               s$N,
                                                               s$corr,
                                                               s$moms,
                                                               s$supp,
                                                               marg.model = j,
                                                               variable.names = s$names,
                                                               checkdata = TRUE,
                                                               tabulate.similar.data = TRUE,
                                                               assume.all.smooth = FALSE,
                                                               s.seed = jiseed
                                                               )
                        )
                    )
                    return(artificial.data)
                }
                )
                names(out) <- trialnames
                return(out)
            }
        )
    )
    names(artds) <- c("gamma", "johnson")
    return(artds)

}
