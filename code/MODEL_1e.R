## MODEL_1e.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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

## model 1.e MHR

model_1e <- c("Death" = fos_1e, "Disease" = fdfs_1e)

#### RUN MODEL

system.time(
RES_model_1e <- lapply(model_1e, function(m)
    lapply(pd, function(d) BootStat(d, "random", reg_formula = m)  ## random = coxme
           )    )
)


## NOTES: 
## median HR increase of ~ 105% for death due to transiting to center with higher risk
## median HR increase of ~ 100% for disease due to transiting to center with higher risk
                              # compare with HR for trt (using model 1.d): does HR lie between 1/MHR and MHR ? If yes: ...
## ... the effect of clustering is bigger than the treatment effect


### PRINT TABLE

t_m_1e <- do.call("rbind", lapply(RES_model_1e, function(m) do.call("cbind", lapply(m, function(d)
    data.frame( x = c(d$Bag$ranefsd1, mhr(d$Bag$ranefsd1), d$Bag$lr, d$Bag$df  )  )   )  )  ) ) 


colnames(t_m_1e) <- names(RES_model_1e[[1]])

t_m_1e <- data.frame("Out." = rep(names(RES_model_1e), rep(dim(t_m_1e)[1]/2, 2)), 
                     "Est." = rep(c("xi", "mhr", "lr", "df"), 2), t_m_1e   )


### TABLE model 1.e FOR SUPPLEMENTARY MATERIAL

t_m_1e_mr <- make.multirow(t_m_1e, 1, rotate = T)

t_m_1e_bold <- cbind(t_m_1e_mr[, 1:2],
                     as.data.frame(t(apply(t_m_1e_mr[, -c(1,2)], 1, function(x)
                         Boldnums(x, tab_thres, abs(x - x[1]), round_up )))) )


print("New Objects: 't_m_1e', 't_m_1e_bold'")






