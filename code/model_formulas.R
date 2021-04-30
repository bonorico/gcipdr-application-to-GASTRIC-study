## model_formulas.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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


##### MODEL FORMULAS for one stage Cox regression

    ## MODEL 1.a (STRATIFIED COX -- EQUALS DISTRIBUTED REGRESSION)

fos_1a <- Surv( timeT, statusT )~trt + strata(trialID)  

fdfs_1a <- Surv( timeS, statusS )~trt + strata(trialID) 

   ## MODEL 1.b (TEST PH ASSUMPTION) 

fos_1b <- Surv( timeT, statusT )~trt + trt:exp(timeT) + strata(trialID)  

fdfs_1b <- Surv( timeS, statusS )~trt + trt:exp(timeS) + strata(trialID) 


    # MODEL 1.c (fixed-effect marginal model)

fos <- Surv( timeT, statusT )~trt               # overall survival 

fdfs <- Surv( timeS, statusS )~trt             # disease free survival surrogate


    # MODEL 1.d (ALL BC VARIABILITY IS CAPTURED BY FRAILTY -- strong assumption) but needed for testing for a center effect 

fos_1d <- Surv( timeT, statusT )~( 1|as.factor(trialID) ) + trt    

fdfs_1d <- Surv( timeS, statusS )~( 1|as.factor(trialID) ) + trt 

    ## MODEL 1.e NULL MODEL (MEDIAN HR)

mhr <- function(logFRsd) exp(logFRsd * qnorm(0.75) * sqrt(2))

fos_1e <- Surv( timeT, statusT )~( 1|as.factor(trialID) )   

fdfs_1e <- Surv( timeS, statusS )~( 1|as.factor(trialID) ) 


    ## MODEL 1.f (THIS IS AKIN TO TWO STAGE MA but in one stage) -- > test for random trt effect 

fos_1f <- Surv( timeT, statusT )~( trt|as.factor(trialID) ) + trt + strata(trialID)

fdfs_1f <- Surv( timeS, statusS )~( trt|as.factor(trialID) ) + trt + strata(trialID)


   ## MODEL 1.g (MULTI LEVEL)

fos_1g <- Surv( timeT, statusT )~( 1 + trt|as.factor(trialID) ) + trt  # overall survival

fdfs_1g <- Surv( timeS, statusS )~( 1 + trt|as.factor(trialID) ) + trt # disease free survival surrogate 

   ## MODEL 1.h (JOINT SURROGATE ENDPOINTS)
