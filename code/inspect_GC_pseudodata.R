## inspect_GC_pseudodata.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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



### extract GC pseudodata results

artdsres <- lapply(artds, function(j) lapply(j, function(x) x$res) )  ## ARTIFICIAL DATA

# read optimization log among centers
artdslog <- lapply(artds, function(j) do.call("rbind", lapply(j, function(x) x$logmex) ) )

 # review optimization procedure (SUPPLEMENT)
table(artdslog$gamma$exitmex)

artdslog$gamma[which(as.character(artdslog$gamma$exitmex) != as.character(artdslog$gamma$bestperf) ), ]
                                        #

table(artdslog$johnson$exitmex)

which(as.character(artdslog$johnson$exitmex) != as.character(artdslog$johnson$bestperf) )

#### check results

lapply(artdsres, function(j) unlist(lapply(j, function(x) class(x)) )) ## center 3 fails optimization GC gamma

## correct center 3 simulation GC gamma

 set.seed( 1608, "L'Ecuyer") 
                 print( system.time(
artdsres[[1]][[3]] <-  with(trial_summaries[["3"]], DataRebuild( H, N, corr, moms, supp, marg.model = "gamma", variable.names = names, checkdata = TRUE, tabulate.similar.data = TRUE)
)
))



# PASTE (POOL) STUDY SPECIFIC ARTIFICIAL DATA (artificial pooled data)

poolartd <-  lapply(artdsres, function(j) PasteDistributed(j, "trialID", trialnames))
names(poolartd) <- c("gamma", "johnson")


print("New object: 'poolartd'")
