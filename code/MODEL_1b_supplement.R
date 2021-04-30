## MODEL_1b_supplement.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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

## ..continuing from MODEL_1b.R

##### SUPPLEMENT 

## inspection of pooled time variable on original IPD and on GC pseudodata

timeT_gc <- lapply(pd[6:7], function(m)
    apply(do.call("rbind", lapply(m, function(d) sort(d$timeT)) ), 2, mean, na.rm = T) ) 


timeS_gc <- lapply(pd[6:7], function(m)
    apply(do.call("rbind", lapply(m, function(d) sort(d$timeS)) ), 2, mean, na.rm = T) ) 


png(file.path(outfolder, "time_transform_model_1b.png"), width = 1000, height = 600)

par(mfrow=c(1, 2)) # approximate shape of timeT variable, motivating "exp" transformation for interaction term

with(gastadj, plot(1:length(timeT), sort(timeT), type = "l",  xlab = "Time rank", ylab = "Time to death (months)", main = "Ranked times to death", col = "red", lwd = 3, ylim = c(0, 300), cex.axis = 1.5, cex.lab = 1.5, cex = 1.5, cex.main = 1.5 ))  

points(1:length(timeT_gc$"GC_few"), sort(timeT_gc$"GC_few"), type = "l", col = "green", lwd = 2)  

points(1:length(timeT_gc$"GC_all"), sort(timeT_gc$"GC_all"), type = "l",  lty = 2, lwd = 2 )  

legend("topleft", c("Orig. IPD", "GC few. mom.", "GC all mom."), lty = c(1,1,2),
       col= c("red", "green", "black"), lwd = c(3,2,2), bty = "n")

with(gastadj, plot(1:length(timeS), sort(timeS), type = "l",  xlab = "Time rank", ylab = "Time to disease (months)", main = "Ranked times to disease", col = "red", lwd = 3, ylim = c(0, 300), cex.axis = 1.5, cex.lab = 1.5, cex = 1.5, cex.main = 1.5 ))  

points(1:length(timeS_gc$"GC_few"), sort(timeS_gc$"GC_few"), type = "l", col = "green", lwd = 2)  

points(1:length(timeS_gc$"GC_all"), sort(timeS_gc$"GC_all"), type = "l",  lty = 2, lwd = 2 )  

dev.off()


#### local IPD inspection for assessment of time-transform for center-specific model 1.b

png(file.path(outfolder, "local_time_transform_DEAT_model_1b.png"), width = 2000, height = 1000)

par(mfrow = c(2, 7))    
for(s in trialnames){
    with(subset(gastadj, subset = trialID == s),
         plot(1:length(timeT), sort(timeT), type = "l", main = s,
              xlab = "Time rank", ylab = "Time to death (months)", col = "red", lwd = 2 ))
}

dev.off()

### disease

png(file.path(outfolder, "local_time_transform_DIS_model_1b.png"), width = 2000, height = 1000)

par(mfrow = c(2, 7))
for(s in trialnames){
    with(subset(gastadj, subset = trialID == s),
         plot(1:length(timeS), sort(timeS), type = "l", main = s,
              xlab = "Time rank", ylab = "Time to disease (months)", col = "red", lwd = 2 ))
}

dev.off()


#### ALTERNATIVE ANALYSES

## NAIVE analysis, accounting for typical privacy-constrained DCN setting, where central computer must a-priori assume a time-transformation equal for each center (e.g. identity), due to lack of IPD detail ...

#### naive distributed regression, assuming identity time-transform for all centers

model_1b_naive <- c("Death" = Surv(timeT, statusT) ~ trt + trt:timeT + strata(trialID),
                    "Disease" = Surv(timeS, statusS) ~ trt + trt:timeS + strata(trialID) )

RES_model_1b_dcn_naive <- lapply(model_1b_naive, function(m)
    lapply(pd[1], function(d) BootStat(d, "fixed", reg_formula = m)  ## fixed = coxph
           )    )

## ...corresponding naive two-stage approach (it wrongly detects a significant non zero interaction ...)..... stage one

ma_model_dcn_naive <- c("Death" = Surv(timeT, statusT) ~ trt + trt:timeT,
                    "Disease" = Surv(timeS, statusS) ~ trt + trt:timeS)

maos_dat_dcn_naive <- lapply(ma_model_dcn_naive, function(m) do.call("rbind", lapply(trialnames, function(s){

 mod <-  stat(subset(gastadj, subset = trialID == s), m, "coxph" )

    out <- data.frame(trialID = s, logHR = mod$coef[2], SDlogHR = mod$sd[2])   } ) ) ) ### NOTE: convergence problems ...

# stage two 
ma_1b_naive <-  lapply(maos_dat_dcn_naive, function(m)
                       metagen(logHR, SDlogHR, trialID, m, sm = "HR", comb.fixed = T))

## print table

t_m_1b_dcn <- na.omit( do.call("rbind", lapply(RES_model_1b_dcn_naive, function(m) do.call("cbind", lapply(m, function(d)
    data.frame( IPD = c(d$Bag$coef, d$Bag$sd, d$Bag$coef/d$Bag$sd, pnorm(d$Bag$coef/d$Bag$sd, lower = T)*2   )  )   )  )  ) ) )

t_m_1b_dcn <- data.frame("Out." = rep(names(RES_model_1b_dcn_naive), rep(dim(t_m_1b_dcn)[1]/2, 2)),
                     "Est." = rep(c("beta", "betaT", "sdbeta", "sdbetaT", "z.1", "z.2", "p.1", "p.2"), 2), t_m_1b_dcn )


t_m_1b_dcn$"fixedMA" <- c(NA, ma_1b_naive$Death$"TE.fixed", NA, ma_1b_naive$Death$"seTE.fixed", NA, ma_1b_naive$Death$"zval.fixed", NA, ma_1b_naive$Death$"pval.fixed", NA, ma_1b_naive$Disease$"TE.fixed", NA, ma_1b_naive$Disease$"seTE.fixed", NA, ma_1b_naive$Disease$"zval.fixed", NA, ma_1b_naive$Disease$"pval.fixed" )

t_m_1b_dcn_mr <- make.multirow(t_m_1b_dcn, 1, rotate = T)

t_m_1b_dcn_bold <- cbind(t_m_1b_dcn_mr[, 1:2],
                         as.data.frame(t(apply(t_m_1b_dcn_mr[, -c(1,2)], 1, function(x)
                             Boldnums(x, tab_thres, abs(x - x[1]), round_up )))) )


## Global (weighted residuals) on MODEL_1a CAN YOU PERFMOR IT IN DCN TODAY ? NO it uses IPD
## FOR DISCUSSION: pval for dfs tends toward 0.05, but this is not detected by CART or GC

model_1a <- c("Death" = fos_1a, "Disease" = fdfs_1a)

system.time(
RES_model_1b_resid <- lapply(model_1a, function(m)
    lapply(pd, function(d) BootStat(d, "residuals", reg_formula = Surv( timeT, statusT )~trt + strata(trialID)  )  ## residuals = cox.zph
           )    )
)

##

t_m_1b_res <- do.call("rbind", lapply(RES_model_1b_resid, function(m) do.call("cbind", lapply(m, function(d) do.call( "rbind", lapply(d$Bag, function(x) x )     ) )  )  ) ) 
           
colnames(t_m_1b_res) <- names(RES_model_1b_resid[[1]])

t_m_1b_res <- data.frame("Out." = c("Death","","Dis.",""), "Est." = rownames(t_m_1b_res), t_m_1b_res   )

t_m_1b_res[t_m_1b_res$"Est." == "pval", -c(1:2)] <- apply(t_m_1b_res[t_m_1b_res$"Est." == "chisq", -c(1:2)], 2, function(x) pchisq(x, 1, lower = F))

t_m_1b_res_bold <- cbind(t_m_1b_res[, 1:2],
                         as.data.frame(t(apply(t_m_1b_res[, -c(1,2)], 1, function(x)
                             Boldnums(x, tab_thres, abs(x - x[1]), round_up )))) )


print("New Objects: t_m_1b_dcn, t_m_1b_dcn_bold, t_m_1b_res, t_m_1b_res_bold")




##   test <- cox.zph(coxph(model_1a[[2]], gastadj, x = TRUE), global = F)
##   test
## plot(test)
##  abline(h = ipdos_1b$coef, lty = 2, col = "red")
## ### FUTURE WORK: graphical residuals analysis on pseudodata
