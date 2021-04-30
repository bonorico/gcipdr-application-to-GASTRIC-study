## MODEL_1h.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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

## model 1.h joint surrogate (two-steps copula)

### work on random subset of original repetitions to speed up calculations

subpd <- lapply(pd[-1], function(d){
    set.seed(673)
    keep <- sample(1:H, 100)   ### just sample 100 reps
    lapply(keep, function(i) d[[i]])    
} ) 

subpd <- c(pd[1], subpd); names(subpd)[1] <- "IPD"


                                        #
getOption("show.error.messages")
options(show.error.messages = FALSE)  ### override internal surrosurv setting in handling errors

system.time(
RES_model_1h <- lapply(subpd, function(d) surroboot(d) )
)

options(show.error.messages = TRUE)


#### PRINT TABLE

t_m_1h <- do.call("cbind", lapply(RES_model_1h, function(d) do.call( "rbind", lapply(d$Bag, function(x) x )     ) )  )   
           
colnames(t_m_1h) <- names(RES_model_1h)

### TABLE model 1.h FOR SUPPLEMENTARY MATERIAL

t_m_1h_bold <- as.data.frame(t(apply(t_m_1h, 1, function(x)
    Boldnums(x, tab_thres, abs(x - x[1]), round_up )))) 

t_m_1h_bold <- data.frame("Est." = rownames(t_m_1h), t_m_1h_bold)

t_m_1h <- data.frame("Est." = rownames(t_m_1h), t_m_1h)

rm(subpd)

print("New Objects: 't_m_1h', 't_m_1h_bold'")



### TODO (discussion): similar to study on cox residuals, future work can also focus on recovery of graphical surrogacy diagnostics ....

## ipdjoint_Cl <- surrosurv(
##     data.frame(id = 1:dim(gastadj)[1], trialref = gastadj[, 1], gastadj[, -1] ),
##     model = "clayton")

## plot(ipdjoint_Cl)


