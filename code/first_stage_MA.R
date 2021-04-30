## first_stage_MA.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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

### first stage of two-stage approach: compute local HRs
### model_1a, model_1f


maos_dat <- do.call("rbind", lapply(trialnames, function(s){

 mod <-  stat(subset(gastadj, subset = trialID == s), fos, "coxph" )

    out <- data.frame(trialID = s, logHR = mod$coef, SDlogHR = mod$sd)


} ) )

                                        #
madfs_dat <- do.call("rbind", lapply(trialnames, function(s){

 mod <-  stat(subset(gastadj, subset = trialID == s), fdfs, "coxph")

    out <- data.frame(trialID = s, logHR = mod$coef, SDlogHR = mod$sd)


} ) )
