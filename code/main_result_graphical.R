## main_results_graphical.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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


                                        # synthesis of all regression results --> are pseudodata-based inferences well correlated with IPD inferences ?

# tab <- readRDS(file.path(outfolder, "all_tabs.rds") )


tmraw <- rbind(
    data.frame(all_tabs$Tab_m1a$raw[, -dim(all_tabs$Tab_m1a$raw)[2]],
               MA = all_tabs$Tab_m1a$raw[, dim(all_tabs$Tab_m1a$raw)[2]],
               Model = "1.a"),
    data.frame(all_tabs$Tab_m1b$raw[, -dim(all_tabs$Tab_m1b$raw)[2]],
               MA = all_tabs$Tab_m1b$raw[, dim(all_tabs$Tab_m1b$raw)[2]],
               Model = "1.b"),
    data.frame(all_tabs$Tab_m1c$raw,
               MA = NA,
               Model = "1.c"),
    data.frame(all_tabs$Tab_m1d$raw,
               MA = NA,
               Model = "1.d"),
    data.frame(all_tabs$Tab_m1e$raw,
               MA = NA,
               Model = "1.e"),
    data.frame(all_tabs$Tab_m1f$raw[, -dim(all_tabs$Tab_m1f$raw)[2]],  
               MA = all_tabs$Tab_m1f$raw[, dim(all_tabs$Tab_m1f$raw)[2]],
               Model = "1.f"),
    data.frame(all_tabs$Tab_m1g$raw,   
               MA = NA,
               Model = "1.g"),
    data.frame(Out. = "Death", all_tabs$Tab_m1h$raw,
               MA = NA,
               Model = "1.h")
)

tmraw <- subset(tmraw, (Est. != "df" & Est. != "conv" & Est. != "coef<mhr"))
tmraw$Est. <- as.factor(tmraw$Est.)
levels(tmraw$Est.)[c(1,2,4,9:11,14:16, 20:21)] <-  c(rep("beta", 3), rep("pval", 3), rep("SDbeta", 3), rep("ztest", 2)) 


### main regression results compact

names_X <- c("Orig. IPD", "IPD boot.", "Dist. boot.", "CART", "VAE", "GC few. mom.", "GC all mom.", "Meta an.")

regr_syn_deat <- data.frame(out = "Death",
                            print_synopsis(
                                subset(tmraw, Out. == "Death", select = 3:10), names_X)
                            )

regr_syn_dis <- data.frame(out = "Disease",
                           print_synopsis(
                               subset(tmraw, Out. == "Disease", select = 3:10), names_X)
                           )

regr_synops <- rbind(regr_syn_deat, regr_syn_dis)

regr_synops_mr <- make.multirow(regr_synops, 1, rotate = T)

print("New Objects: 'regr_synops', 'regr_synops_mr'.")





## # NOT RUN

## tmdat <- tmraw %>% tibble() %>% 
##     pivot_longer(cols = 4:10,
##                  names_to = "Type",
##                  values_to = "approx",
##                  values_ptypes =
##                      list(approx = numeric())
##                  )




## tmdat %>% subset(Est. == "Chi", select = IPD) %>% summary()
## tmdat %>% subset(Est. == "lr", select = IPD) %>% summary()

##                                         # resize Chi and LR values grater than 10 for comparability
## subs <- which(tmdat$Est. == "Chi" | tmdat$Est. == "lr")
## reg_dat <- tmdat

## reg_dat[subs, ] <- reg_dat[subs, ] %>%
##     mutate(.keep = "unused",
##                     IPD = if_else(IPD>10, log(IPD), IPD),
##                     approx = if_else(approx>10, log(approx), approx))


## hist(reg_dat$IPD)

## apprx_GC_few_deat <- lm(IPD~approx, data = reg_dat,
##                         subset = (Out. == "Death" & Type == "GC_few")
##                         )

## apprx_GC_all_deat <- lm(IPD~approx, data = reg_dat,
##                         subset = (Out. == "Death" & Type == "GC_all")
##                         )



## summary(apprx_GC_few_deat)
## summary(apprx_GC_all_deat)

## par(mfrow=c(1, 2))
## plot(IPD~approx, data = reg_dat, main = "few. mom.",
##      subset = (Out. == "Death" & Type == "GC_few"))
## abline(apprx_GC_few_deat, col = "red")
## plot(IPD~approx, data = reg_dat, main = "all mom.",
##      subset = (Out. == "Death" & Type == "GC_all"))
## abline(apprx_GC_all_deat, col = "red")








