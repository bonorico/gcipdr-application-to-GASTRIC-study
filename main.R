## main.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository).
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


source("load_dependencies.R")

options(mc.cores = 3L)   # to increase parallel cores

## to locate warning: options(warn = 2) ## default warn = 0

outfolder <- "./results"

## BLOCK 1 #### load GASTRIC IPD

if( !any( as.numeric(table(gastadj$patientID)) > 1) ){
    gastadj <- gastadj[, -2]

  }else{
    print("There are patients with multple obss. Keeping var 'patientID'")
  }

  gastadj[, c(3,5)] <-  (gastadj[, c(3,5)]/365)*12  # monthly time scale CORRECT

  trialnames <- levels(as.factor(gastadj$trialID))


## BLOCK 2 ##### EXTRACT key IPD summaries

trial_summaries <- lapply(trialnames, function(s){

  dat <- subset(gastadj, subset = trialID == s, select = -1)

    datin <- Return.key.IPD.summaries( Return.IPD.design.matrix(dat) )

    return(list( moms = datin$first.four.moments, corr = datin$correlation.matrix,
      N = datin$sample.size, supp = datin$is.binary.variable, names = datin$variable.names ))

       }  ); names(trial_summaries) <- trialnames


## BLOCK 3 ##### PSEUDODATA #####################

H <- 300  # set simulation size

## BLOCK 4 ###### GENERATE GC PSEUDODATA

source("code/generate_GC_pseudodata.R")

system.time(
artds <- GC_pseudodata()
)

## BLOCK 5 ###### CHECK AND CORRECT GC PSEUDODATA

source("code/inspect_GC_pseudodata.R", echo = T)

## BLOCK 6 ######## COMPARATIVE PSEUDODATA METHODS (both disclosive and not)

source("code/other_pseudodata_methods.R", echo = T)

### all data in a list (convenient to run models in loop)

pd <- list("IPD" = list(gastadj[, c(2:6, 1)]),
           "IPDboot" = lapply(bootres, function(x) x[, c(2:6, 1)]),
           "Distboot" = bootresDist,
           "CART" = pool.synpdat,
           "VAE" = pool.vaedat,
           "GC_few" = poolartd$gamma,
           "GC_all" = poolartd$johnson
           )

rm(bootres, bootresDist, pool.synpdat, pool.vaedat)  # a bit clean-up

# saveRDS(pd, "all_pseudodata.rds")
#  pd <- readRDS(file.path(outfolder, "all_pseudodata.rds"))

## BLOCK 7 ##### MODEL FORMULAS

source("code/model_formulas.R")

## BLOCK 8 ### for TWO-STAGE approach
source("code/first_stage_MA.R", echo = T)

## BLOCK 9 ##### RUN ONE-STAGE REGRESSION ANALYSES ON DATA #####

round_up <- 3
tab_thres <- 0.05  ## 0.01

source("code/MODEL_1a.R", echo = T)

source("code/MODEL_1b.R", echo = T)

source("code/MODEL_1b_supplement.R", echo = T)

source("code/MODEL_1c.R", echo = T)

source("code/MODEL_1d.R", echo = T)

source("code/MODEL_1e.R", echo = T)

source("code/MODEL_1f.R", echo = T)

source("code/MODEL_1g.R", echo = T)

source("code/MODEL_1h.R", echo = T)

## BLOCK 10 #### ORIGINAL IPD KM ANALYSES

source("code/IPD_KM.R", echo = T)
source("code/IPD_KM_plots.R", echo = T)  ## plot curves

## BLOCK 11 #### CENTER-SPECIFIC BASELINE HAZARDS (CENTER EFFECT)

source("code/compute_center_KMs_trt.R", echo = T)  ### center effect on KM under trt
source("code/compute_center_BHs.R", echo = T)  ### BH proportionality

pdf(file.path(outfolder, "KM_trt_BH_centers_plot_death_John.eps"), width = 10.4, height = 12.4)
plot_grid(main3, main4, labels="AUTO", nrow = 2, rel_heights = c(1.2, 1))
dev.off()
## BLOCK 12 ######## SUBGROUP ANALYSIS #########

source("code/subgroup_analysis.R", echo = T)

## BLOCK 13 ######## DESCRIPTIVE STATISTICS

source("code/pseudodata_descriptive.R", echo = T)

## BLOCK 14 #### MORE DATA CHECKING FOR APPENDIX

source("code/appendix_res_1.R", echo = T)

#### TOREPORT all results

all_tabs <- list(
     "Tab_m1a" = list("raw" = t_m_1a, "latex" = t_m_1a_bold),
     "Tab_m1b" = list("raw" = t_m_1b, "latex" = t_m_1b_bold),
     "Tab_m1c" = list("raw" = t_m_1c, "latex" = t_m_1c_bold),
     "Tab_m1d" = list("raw" = t_m_1d, "latex" = t_m_1d_bold),
     "Tab_m1e" = list("raw" = t_m_1e, "latex" = t_m_1e_bold),
     "Tab_m1f" = list("raw" = t_m_1f, "latex" = t_m_1f_bold),
     "Tab_m1g" = list("raw" = t_m_1g, "latex" = t_m_1g_bold),
     "Tab_m1h" = list("raw" = t_m_1h, "latex" = t_m_1h_bold),
    "Tab_m1b_dcn" = list("raw" = t_m_1b_dcn, "latex" = t_m_1b_dcn_bold),  ## distributed approach
    "Tab_m1b_res" = list("raw" = t_m_1b_res, "latex" = t_m_1b_res_bold),  ## Schönfeld residuals test
    "Tab_2_subg" = list("raw" = table2, "latex" = table2_bold),
    gastadj = gastadj, descr1 = descr1, descr2 = descr2,
    mom_dat = mom_dat_bold, corr_dat = corr_dat_bold,
    maos_1f = maos_1f, madfs_1f = madfs_1f
)


source("code/main_result_graphical.R")  # regression results synopsis (a high level synthesis)

all_tabs$regr_synops <- list("raw"= regr_synops, "latex" = regr_synops_mr )
all_tabs$km_synops <- list("raw"= km_synops_out, "latex" = km_synops_out_mr )

saveRDS(all_tabs, file.path(outfolder, "all_tabs.rds"))



Sweave("results/main_doc.Rnw", encoding="UTF8"); tools::texi2pdf("results/main_doc.tex")
Sweave("results/supplementary_material.Rnw", encoding="UTF8"); tools::texi2pdf("results/supplementary_material.tex")
