## IPD_KM.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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


##### TWO-STAGE SURV (parametric survival) based on MLE haz
                                        # note can tweak kmparam to 'random.effect = TRUE' to obtain the alternative estimate

km.MAos <- data.frame(Outcome = "Overall survival probability", Estimate = "Meta an.",
                      with(gastadj, kmparam(timeT, statusT, trt, 300, trialID)) )
km.MAdfs <- data.frame(Outcome = "Disease free survival probability", Estimate = "Meta an.",
                       with(gastadj, kmparam(timeS, statusS, trt, 300, trialID)) )


#### ONE-STAGE APPROACHES

model_km <- c("Death" = fos, "Disease" = fdfs)

#### RUN MODEL

 system.time(
RES_model_km <- lapply(model_km, function(m)
    lapply(pd, function(d) BootStat(d, "km", reg_formula = m)  ## fixed = coxph
           )    )
)

#
kmdata <- do.call("rbind", lapply(RES_model_km, function(m) do.call("rbind", lapply(m, function(d) d$Bag  )) ) )

kmdata <- data.frame(
    Outcome = as.factor(ifelse(
        unlist(lapply(rownames(kmdata), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1] ) ) == "Death",
        "Overall survival probability", "Disease free survival probability") ),
    Estimate = as.factor(unlist(lapply(rownames(kmdata), function(x) strsplit(x, ".", fixed = TRUE)[[1]][2] ) ) ),
    kmdata
        )

kmdata <- rbind(kmdata, km.MAos, km.MAdfs)

levels(kmdata$Estimate)[2:6] <- c("Dist. boot.", "GC all mom.", "GC few. mom.", "Orig. IPD", "IPD boot.")
kmdata$Estimate <- factor(kmdata$Estimate, levels = c("Orig. IPD", "IPD boot.", "Dist. boot.", "CART",
                                                      "VAE", "GC few. mom.", "GC all mom.", "Meta an."))

kmdata$trt <- ifelse(kmdata$trt == 0, "Control", "Treatment")

colnames(kmdata)[7] <- "Treatment"

kmdata$Outcome <- relevel(kmdata$Outcome, "Overall survival probability")


### main summary --> Kolmogorov-Smirnov test (unweighted for censoring)

km_synops <- do.call("rbind",
                     lapply(levels(kmdata$Outcome),
                            function(o)
                                data.frame(
                                    out = o,
                                    do.call("cbind",
                                            lapply(unique(kmdata$Treatment),
                                                   function(g)
                                                       as.data.frame(                                                       
                                                           km_ks_test(
                                                               subset(kmdata,
                                                                      Outcome == o &
                                                                      Treatment == g
                                                                      )
                                                               
                                                           )
                                                       )
                                                   
                                                   )
                                            )[, -5]   # drop double var
                                )
                            )
                     )

    
km_synops$out <- as.factor(km_synops$out)
levels(km_synops$out) <- c("DFS", "OS")
names(km_synops)[3:5] <- paste0("CRT_", names(km_synops)[3:5])
names(km_synops)[6:8] <- paste0("TRT_", names(km_synops)[6:8])

km_synops_compact <- (km_synops[, 3:4] + km_synops[, 6:7])/2

km_synops_out <- data.frame(
    km_synops[, 1:2],
    km_synops_compact,
    as.integer(km_synops[, 5] & km_synops[, 8])
)
names(km_synops_out)[-c(1:2)] <- c("KS_stat", "sig_lev", "H0_refused")

all( (km_synops_out[, 3] > km_synops_out[, 4]) == km_synops_out[, 5])  # check all TRUE

km_synops_out_mr <- make.multirow(km_synops_out, 1, rotate = T)


print("New Objects: 'km_synops_out', 'km_synops_out_mr'.")

# km_synops[, 3:8] <- round(km_synops[, 3:8], 3)











### TO DELETE
## # NOTE: Rsq can be negative for KM values since KM doeas not obviously have to be linear
## km_synops <- lapply(kmwide, function(x)
##     lapply(x,
##            function(y)
##                lapply(y,
##                       function(j)
##                           print_synopsis(j)
##                       )
##            )
##     )

##                                         # overall synopsis is an average of synopsis for y and x axis and over trt
## syn_names <- km_synops[[1]][["trt"]][["time"]]$pair
## km_synops_aver <- lapply(km_synops,
##                          function(x)
##                              lapply(trgt[-2],
##                                     function(j)
##                                         data.frame(
##                                             pair = syn_names,
##                                                                        # double matrix sum over time and treatment group
                                                
##                                              (       ((
##                                                         x[["ctr"]][[j]][-1] +
##                                                         x[["ctr"]][["time"]][-1] )/2)
##                                             +
##                                                     ((
##                                                         x[["trt"]][[j]][-1] +
##                                                         x[["trt"]][["time"]][-1] )/2) 
##                                             )/2
##                                         )

##                                     )
##                          )


## for (o in names(km_synops_aver))
##     names(km_synops_aver[[o]]) <- trgt[-2]


## data.frame(km_synops_aver[[1]][["km"]][1],
##            round(km_synops_aver[[1]][["km"]][-1], 3)
##            )[, -3]

## plot(kmwide[[1]][["trt"]][["time"]]$"IPD boot",
##      kmwide[[1]][["trt"]][["time"]]$"Orig. IPD")

