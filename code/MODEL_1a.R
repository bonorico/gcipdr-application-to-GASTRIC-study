## MODEL_1a.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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


## model 1.a (stratification by center)

model_1a <- c("Death" = fos_1a, "Disease" = fdfs_1a)

#### RUN MODEL

system.time(
RES_model_1a <- lapply(model_1a, function(m)
    lapply(pd, function(d) BootStat(d, "fixed", reg_formula = m)  ## fixed = coxph
           )    )
)


#### TWO STAGE FIXED EFFECTS

# second stage                                        
                                        # # fixed-effect meta-analysis

maos_1a <-  metagen(logHR, SDlogHR, trialID, maos_dat, sm = "HR", comb.fixed = T)
#
madfs_1a <-  metagen(logHR, SDlogHR, trialID, madfs_dat, sm = "HR", comb.fixed = T)


### PRINT TABLE


t_m_1a <- na.omit( do.call("rbind", lapply(RES_model_1a, function(m) do.call("cbind", lapply(m, function(d) do.call( "rbind", lapply(d$Bag, function(x) x )     ) )  )  ) ) )
           
colnames(t_m_1a) <- names(RES_model_1a[[1]])

t_m_1a <- data.frame("Out." = rep(names(RES_model_1a), rep(dim(t_m_1a)[1]/2, 2)), "Est." = rownames(t_m_1a), t_m_1a   )

t_m_1a$"fixedMA" <- c(maos_1a$TE.fixed, maos_1a$seTE.fixed, NA, NA, madfs_1a$TE.fixed, madfs_1a$seTE.fixed, NA, NA)

### add beta z-test

t_m_1a <- ZtestL(t_m_1a)

### TABLE model 1.a FOR SUPPLEMENTARY MATERIAL

t_m_1a_mr <- make.multirow(t_m_1a, 1, rotate = T)

t_m_1a_bold <- cbind(t_m_1a_mr[, 1:2],
                     as.data.frame(t(apply(t_m_1a_mr[, -c(1,2)], 1, function(x)
                         Boldnums(x, tab_thres, abs(x - x[1]), round_up) ) ) ) )


print("New Objects: 't_m_1a', 't_m_1a_bold'")







    


 

   
    
