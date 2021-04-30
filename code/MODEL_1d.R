## MODEL_1d.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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


## model 1.d TEST FOR CENTER EFFECT

model_1d <- c("Death" = fos_1d, "Disease" = fdfs_1d)

#### RUN MODEL

system.time(
RES_model_1d <- lapply(model_1d, function(m)
    lapply(pd, function(d) BootStat(d, "random", reg_formula = m)  ## random = coxme
           )    )
)



### PRINT TABLE

t_m_1d <- do.call("rbind", lapply(RES_model_1d, function(m) do.call("cbind", lapply(m, function(d)
    data.frame( x = c(d$Bag$coef, d$Bag$sd, d$Bag$ranefsd1,
                      mhr(d$Bag$ranefsd1),
                      (exp(d$Bag$coef) > 1/mhr(d$Bag$ranefsd1)) & (exp(d$Bag$coef) < mhr(d$Bag$ranefsd1)),
                      d$Bag$lr, d$Bag$df  )  )   )  )  ) ) 



colnames(t_m_1d) <- names(RES_model_1d[[1]])

t_m_1d <- data.frame("Out." = rep(names(RES_model_1d), rep(dim(t_m_1d)[1]/2, 2)),
                     "Est." = rep(c("coef", "sd", "xi", "mhr", "coef<mhr", "lr", "df"), 2), t_m_1d   )

lr_d <- subset(t_m_1d, Est. == "lr" | Est. == "df")
lr_c <- subset(t_m_1c, Est. == "lr" | Est. == "df")

### LR model_d vs model_c (test center effect)

lr <- lapply(c("Death", "Disease"), function(o)
    LRtest(as.numeric(subset(lr_d, Out. == o & Est. == "lr", select = -c(1:2))),
           as.numeric(subset(lr_c, Out. == o & Est. == "lr", select = -c(1:2))),
           as.numeric(subset(lr_d, Out. == o & Est. == "df", select = -c(1:2))),
           as.numeric(subset(lr_c, Out. == o & Est. == "df", select = -c(1:2)) ))
       )

lr_death <- as.data.frame(t(lr[[1]]))  
lr_dis <- as.data.frame(t(lr[[2]]))
colnames(lr_death) <- colnames(lr_dis) <- names(pd)

t_m_1d <- rbind(subset(t_m_1d, Out. == "Death"),
                data.frame(Out. = "Death", Est. = rownames(lr_death), lr_death),
                subset(t_m_1d, Out. == "Disease"),
                data.frame(Out. = "Disease", Est. = rownames(lr_dis), lr_dis) )
                
## add beta z-test

t_m_1d <- ZtestL(t_m_1d)

### TABLE model 1.d FOR SUPPLEMENTARY MATERIAL

t_m_1d_mr <- make.multirow(t_m_1d, 1, rotate = T)

t_m_1d_bold <- cbind(t_m_1d_mr[, 1:2],
                     as.data.frame(t(apply(t_m_1d_mr[, -c(1,2)], 1, function(x)
                         Boldnums(x, tab_thres, abs(x - x[1]), round_up )))) )



print("New Objects: 't_m_1d', 't_m_1d_bold'")





