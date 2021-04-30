## pseudodata_descriptive.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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


######## DESCRIPTIVE STATISTICS

descr1 <- summ_summs(trial_summaries, T, rd = 2)

descr2 <- summ_summs(trial_summaries, T, T, rd = 2)



## ## GOODNESS OF ARTIFICIAL DATA
### tables 

artsumm <- lapply(
    pd
    ## list("Orig. IPD" = list(
    ##          gastadj[, c(3:7, 1)]
    ##      ),
    ##      "IPD boot." = lapply(bootres,
    ##                           function(x)
    ##                               x[, c(3:7, 1)]
    ##                           ),
    ##      "Dist. boot." = lapply(bootresDist,
    ##                           function(x)
    ##                               x[, -1]
    ##                          ),
    ##      "CART" = pool.synpdat,
    ##      "VAE" = pool.vaedat,
    ##      "GC few. mom." = poolartd$gamma,
    ##      "GC all mom." = poolartd$johnson
    ##      )
   ,
    function(d)
        diagnostic_summ(d, 6)
) 


mom_dat <- do.call(
    "rbind",
    lapply(names(artsumm[[1]][-5]),
           function(s)
               do.call(
                   "rbind",
                   lapply(names(artsumm),
                          function(d)
                              data.frame(
                                  stat = s,
                                  data = d,
                                  data.frame(
                                      t( artsumm[[d]][[s]])
                                  )
                              )
                          )
               )
           )
)

corr_dat <-  do.call(
    "cbind",
    lapply(artsumm,
           function(d)
               data.frame( d[["corr"]])
           )
)

colnames(corr_dat)  <- names(artsumm)
       

corr_dat_bold <- as.data.frame(
    t(
        apply(corr_dat,
              1,
              function(x)
                  Boldnums(x,
                           0.1,
                           abs(x - x[1])
                           )
              )
    )
)
colnames(corr_dat_bold) <- colnames(corr_dat)


treshs <- c(NA, NA, 0.1, 10, 0.1, 10, 0.1)

mom_dat_bold <- do.call(
    "rbind",
    lapply(levels(as.factor(mom_dat$stat)),
           function(s)
               cbind(mom_dat[mom_dat$stat == s, 1:2],
                     do.call(
                         "cbind",
                         lapply(3:7,
                                function(v)
                                {
                                    y <-  mom_dat[mom_dat$stat == s, v]
                                    Boldnums(y,
                                             ifelse(
                                                 any(s == c("skx", "kdx")),
                                                 treshs[v],
                                                 0.5),
                                             abs(y-y[1])
                                             )

                                }
                                )
                     )
                     )
           )
)

colnames(mom_dat_bold) <- colnames(mom_dat)


levels(mom_dat_bold$stat) <- c("Mean", "Sqrt. Variance", "Skewness", "Kurtosis")

mom_dat_bold <- make.multirow(mom_dat_bold, 1, rotate = T)


print("New Objects: 'descr1', 'descr2', 'mom_dat_bold' 'corr_dat_bold'")
