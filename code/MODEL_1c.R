## MODEL_1c.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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


 ## model 1.c marginal effect

model_1c <- c("Death" = fos, "Disease" = fdfs)

#### RUN MODEL

system.time(
RES_model_1c <- lapply(model_1c, function(m)
    lapply(pd, function(d) BootStat(d, "fixed", reg_formula = m)  ## fixed = coxph
           )    )
)



### PRINT TABLE


t_m_1c <- na.omit( do.call("rbind", lapply(RES_model_1c, function(m) do.call("cbind", lapply(m, function(d) do.call( "rbind", lapply(d$Bag, function(x) x )     ) )  )  ) ) )
           
colnames(t_m_1c) <- names(RES_model_1c[[1]])

t_m_1c <- data.frame("Out." = rep(names(RES_model_1c), rep(dim(t_m_1c)[1]/2, 2)), "Est." = rownames(t_m_1c), t_m_1c   )

## add beta z-test

t_m_1c <- ZtestL(t_m_1c)

### TABLE model 1.c FOR SUPPLEMENTARY MATERIAL

t_m_1c_mr <- make.multirow(t_m_1c, 1, rotate = T)

t_m_1c_bold <- cbind(t_m_1c_mr[, 1:2],
                     as.data.frame(t(apply(t_m_1c_mr[, -c(1,2)], 1, function(x)
                         Boldnums(x, tab_thres, abs(x - x[1]), round_up )))) )


print("New Objects: 't_m_1c', 't_m_1c_bold'")




