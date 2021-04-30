## MODEL_1g.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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

## model 1.g multi-level

model_1g <- c("Death" = fos_1g, "Disease" = fdfs_1g)

#### RUN MODEL

system.time(
RES_model_1g <- lapply(model_1g, function(m)
    lapply(pd, function(d) BootStat(d, "random", reg_formula = m)  ## random = coxme
           )    )
)


## notes:
## mhr on tau higher in multi-levels  (see deJOng for interpreation of MHR here)
## very high interrelatoin between center and treatment-by-center effects !!!


### PRINT TABLE


t_m_1g <-  do.call("rbind", lapply(RES_model_1g, function(m) do.call("cbind", lapply(m, function(d)
    data.frame( x = c(d$Bag$coef, d$Bag$sd, d$Bag$ranefsd1, d$Bag$ranefsd2,
                      d$Bag$revcov, d$BagSE$revcov, d$Bag$lr, d$Bag$df  )  )   )  )  ) ) 

           
colnames(t_m_1g) <- names(RES_model_1g[[1]])

t_m_1g <- data.frame("Out." = rep(names(RES_model_1g), rep(dim(t_m_1g)[1]/2, 2)),
                     "Est." = c("coef", "sd", "xi", "tau", "rho", "mcMSErho", "lr", "df"), t_m_1g   )


lr_g <- subset(t_m_1g, Est. == "lr" | Est. == "df")
lr_f <- subset(t_m_1f, Est. == "lr" | Est. == "df")[, -10]

### LR model_g vs model_f (test center effect)

lr <- lapply(c("Death", "Disease"), function(o)
    LRtest(as.numeric(subset(lr_g, Out. == o & Est. == "lr", select = -c(1:2))),
           as.numeric(subset(lr_f, Out. == o & Est. == "lr", select = -c(1:2))),
           as.numeric(subset(lr_g, Out. == o & Est. == "df", select = -c(1:2))),
           as.numeric(subset(lr_f, Out. == o & Est. == "df", select = -c(1:2))[1,] ))  ## first df    
       )


lr_death <- as.data.frame(t(lr[[1]]))  
lr_dis <- as.data.frame(t(lr[[2]]))
colnames(lr_death) <- colnames(lr_dis) <- names(pd)

t_m_1g <- rbind(subset(t_m_1g, Out. == "Death"),
                data.frame(Out. = "Death", Est. = rownames(lr_death), lr_death),
                subset(t_m_1g, Out. == "Disease"),
                data.frame(Out. = "Disease", Est. = rownames(lr_dis), lr_dis) )

## add beta z-test

t_m_1g <- ZtestL(t_m_1g)

### TABLE model 1.g FOR SUPPLEMENTARY MATERIAL

t_m_1g_mr <- make.multirow(t_m_1g, 1, rotate = T)

t_m_1g_bold <- cbind(t_m_1g_mr[, 1:2],
                     as.data.frame(t(apply(t_m_1g_mr[, -c(1,2)], 1, function(x)
                         Boldnums(x, tab_thres, abs(x - x[1]), round_up )))) )

t_m_1g_bold[, -c(1:3)] <- apply(t_m_1g_bold[, -c(1:3)], 2, as.character)  ## add MSE

t_m_1g_bold[c(8, 21), -c(1:3)] <- t(apply(t_m_1g[c(8, 21), -c(1:3)], 1, function(x) as.character(round(x, 3))))

t_m_1g_bold[, -c(1:3)] <- apply(t_m_1g_bold[, -c(1:3)], 2, as.factor)

print("New Objects: 't_m_1g', 't_m_1g_bold'")




