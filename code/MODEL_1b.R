## MODEL_1b.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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

### model 1.b


#### NON-NAIVE ANALYSIS accounting for empirical inspection of IPD time variable, HAD WE HAD IPD ... (see NAIVE approach in SUPP)

model_1b <- c("Death" = fos_1b, "Disease" = fdfs_1b)

#### RUN MODEL with exp() transformation, according to visual inspection of pseudodata (see SUPPLEMENT)

system.time(
RES_model_1b <- lapply(model_1b, function(m)
    lapply(pd, function(d) BootStat(d, "fixed", reg_formula = m)  ## fixed = coxph
           )    )
)


### TWO-STAGE APPROACH to estimate interaction effect: each center inspects IPD locally (see SUPPLEMENT below) and custom-model time interaction (it gives correct but little reliable results, because local models struggle to converge throwing warnings)

### stage one

 ff_os <- paste0("Surv(timeT, statusT) ~ trt + trt:",
                   c("exp(timeT)", rep("log(timeT)",3), rep("exp(timeT)", 2), "timeT",
                     "exp(timeT)", "timeT", rep("exp(timeT)", 3), "timeT", "exp(timeT)") )  

  names(ff_os) <- trialnames

maos_dat_ph_non_naive <- do.call("rbind", lapply(trialnames, function(s){

 mod <-  stat(subset(gastadj, subset = trialID == s), as.formula(ff_os[s]), "coxph" )

    out <- data.frame(trialID = s, logHR = mod$coef[2], SDlogHR = mod$sd[2])   } ) )  ### NOTE: convergence problems ...

                                        #  disease 
                                        
 ff_dfs <- paste0("Surv(timeS, statusS) ~ trt + trt:",
                  c("exp(timeS)", rep("log(timeS)", 3), rep("exp(timeS)", 2),
                    "timeS", "exp(timeS)", "timeS", rep("exp(timeS)", 3), "timeS", "exp(timeS)"  ) )  

  names(ff_dfs) <- trialnames

                                        #

madfs_dat_ph_non_naive <- do.call("rbind", lapply(trialnames, function(s){

 mod <-  stat(subset(gastadj, subset = trialID == s),  as.formula(ff_dfs[s]), "coxph")

    out <- data.frame(trialID = s, logHR = mod$coef[2], SDlogHR = mod$sd[2])   } ) )

### stage two: Ratio of largest to smallest sampling variance extremely large. May not be able to obtain stable results.

maos_1b_non_naive <-  metagen(logHR, SDlogHR, trialID, maos_dat_ph_non_naive, sm = "HR", comb.fixed = T)
#
madfs_1b_non_naive <-  metagen(logHR, SDlogHR, trialID, madfs_dat_ph_non_naive, sm = "HR", comb.fixed = T)


### PRINT TABLE

t_m_1b <- do.call("rbind", lapply(RES_model_1b, function(m) do.call("cbind", lapply(m, function(d)
    data.frame( x = c(d$Bag$coef, d$Bag$sd, d$Bag$coef/d$Bag$sd, pnorm(d$Bag$coef/d$Bag$sd, lower = T)*2   )  )   )  )  ) ) 
           
colnames(t_m_1b) <- names(RES_model_1b[[1]])

t_m_1b <- data.frame("Out." = rep(names(RES_model_1b), rep(dim(t_m_1b)[1]/2, 2)),
                     "Est." = rep(c("beta", "betaT", "sdbeta", "sdbetaT", "z.1", "z.2", "p.1", "p.2"), 2), t_m_1b   )


t_m_1b$"fixedMA" <- c(NA, maos_1b_non_naive$"TE.fixed", NA, maos_1b_non_naive$"seTE.fixed", NA, maos_1b_non_naive$"zval.fixed", NA, maos_1b_non_naive$"pval.fixed", NA, madfs_1b_non_naive$"TE.fixed", NA, madfs_1b_non_naive$"seTE.fixed", NA, madfs_1b_non_naive$"zval.fixed", NA, madfs_1b_non_naive$"pval.fixed" )


t_m_1b_mr <- make.multirow(t_m_1b, 1, rotate = T)

t_m_1b_bold <- cbind(t_m_1b_mr[, 1:2],
                     as.data.frame(t(apply(t_m_1b_mr[, -c(1,2)], 1, function(x)
                         Boldnums(x, tab_thres, abs(x - x[1]), round_up )))) )


print("New Objects: t_m_1b, t_m_1b_bold")

### to be continued in supplement ...



