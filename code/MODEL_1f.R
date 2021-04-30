## MODEL_1f.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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


## model 1.f (trt random effect, aking to 2stage MA)

model_1f <- c("Death" = fos_1f, "Disease" = fdfs_1f)

#### RUN MODEL

system.time(
RES_model_1f <- lapply(model_1f, function(m)
    lapply(pd, function(d) BootStat(d, "random", reg_formula = m)  ## random = coxme
           )    )
)


## alternative method: classic MA (e.g., TWO-STAGE MA)

                                        # # Sidik-Jonkman estimator can best detect variability between centers when true tau is high .... no HK correction used here

maos_1f <-  metagen(logHR, SDlogHR, trialID, maos_dat, sm = "HR", method.tau = "SJ", comb.random = T)
#
madfs_1f <-  metagen(logHR, SDlogHR, trialID, madfs_dat, sm = "HR", method.tau = "SJ", comb.random = T)
                                        #

### MA total variance computation (DER SIMONIAN & LAIRD)

w <- 1/(maos_1f$seTE^2 + maos_1f$tau2)  ### weight
musd <- 1/sqrt(sum(w))   ## seTE.ransom

musd; maos_1f$seTE.random  ## OK

exp(maos_1f$TE.random - qnorm(0.025, lower = F)*maos_1f$seTE.random)
exp(maos_1f$TE.random + qnorm(0.025, lower = F)*maos_1f$seTE.random)

### NAIVE assessment of HR heterogeneity: simple SD
sd(maos_dat$logHR)
sd(madfs_dat$logHR)

##
meta::forest(maos_1f)
f1 <- grid.grab()

meta::forest(madfs_1f) 
f2 <- grid.grab()

png(file.path(outfolder, "forestplots.png"), width = 630, height = 730)

grid.arrange(f1, f2, ncol=1)

dev.off()


### PRINT TABLE


t_m_1f <- do.call("rbind", lapply(RES_model_1f, function(m) do.call("cbind", lapply(m, function(d)
    data.frame( x = c(d$Bag$coef, d$Bag$sd, d$Bag$ranefsd1, d$Bag$lr, d$Bag$df  )  )   )  )  ) ) 

           
colnames(t_m_1f) <- names(RES_model_1f[[1]])

t_m_1f <- data.frame("Out." = rep(names(RES_model_1f), rep(dim(t_m_1f)[1]/2, 2)),
                     "Est." = c("coef", "sd", "tau", "lr", "df"), t_m_1f   )

t_m_1f$"randMA" <- NA  ### insert randMA res soon

lr_f <- subset(t_m_1f, Est. == "lr" | Est. == "df")
lr_a <- subset(t_m_1a, Est. == "lr" | Est. == "df")

### LR model_f vs model_a (test center effect)

lr <- lapply(c("Death", "Disease"), function(o)
    LRtest(as.numeric(subset(lr_f, Out. == o & Est. == "lr", select = -c(1:2))),
           as.numeric(subset(lr_a, Out. == o & Est. == "lr", select = -c(1:2))),
           as.numeric(subset(lr_f, Out. == o & Est. == "df", select = -c(1:2))),
           as.numeric(subset(lr_a, Out. == o & Est. == "df", select = -c(1:2)) ))
       )

lr_death <- as.data.frame(t(lr[[1]]))  
lr_dis <- as.data.frame(t(lr[[2]]))
colnames(lr_death) <- colnames(lr_dis) <- c(names(pd), "randMA")

t_m_1f <- rbind(subset(t_m_1f, Out. == "Death"),
                data.frame(Out. = "Death", Est. = rownames(lr_death), lr_death),
                subset(t_m_1f, Out. == "Disease"),
                data.frame(Out. = "Disease", Est. = rownames(lr_dis), lr_dis) )


t_m_1f$"randMA" <- c(maos_1f$TE.random, maos_1f$seTE.random, maos_1f$tau, rep(NA, 4), maos_1f$pval.Q,
                     madfs_1f$TE.random, madfs_1f$seTE.random, madfs_1f$tau, rep(NA, 4), madfs_1f$pval.Q)

## add beta z-test

t_m_1f <- ZtestL(t_m_1f)

### TABLE model 1.f FOR SUPPLEMENTARY MATERIAL

t_m_1f_mr <- make.multirow(t_m_1f, 1, rotate = T)

t_m_1f_bold <- cbind(t_m_1f_mr[, 1:2],
                     as.data.frame(t(apply(t_m_1f_mr[, -c(1,2)], 1, function(x)
                         Boldnums(x, tab_thres, abs(x - x[1]), round_up  )))) )

## im MA test for tau uses standard Heterogeneity Q-statistic (not comparable with IPD test)
t_m_1f_bold$"randMA" <- as.character(t_m_1f_bold$"randMA")
t_m_1f_bold$"randMA"[c(8:9,18:19)] <- as.character(round( c(maos_1f$Q, maos_1f$df.Q, madfs_1f$Q, madfs_1f$df.Q), 2) )
t_m_1f_bold$"randMA" <- as.factor(t_m_1f_bold$"randMA")


print("New Objects: 't_m_1f', 't_m_1f_bold'")



