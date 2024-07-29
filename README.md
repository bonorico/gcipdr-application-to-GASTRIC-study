# gcipdr-application-to-GASTRIC-study

Run file "main.R" for main analysis.

Run files "sim_study*.R" for simulation analyses (beware each file can take up to 24 hours to execute on 2 cores). Simulations were made with R 3.6.3 under a Linux 5.4.0-81-generic 18.04.1-Ubuntu system on a Intel(R) Core(TM) i7-2640M CPU @ 2.80GHz machine.  

This analysis was presented in [Survival analysis without sharing of individual patient data by using a Gaussian copula](https://onlinelibrary.wiley.com/doi/10.1002/pst.2415)

# Detail

Tensorflow installation for VAE might be problematic. The code attempts to distinguish between OS (Windows or not) and tailor installation accordingly. Because of rapid development of packages tensorflow and keras the current VAE code might become outdated any time soon. Efforts to mantain the VAE code cannot be guaranteed.
VAE analyses could be lately run on python 3.10.11 (via reticulate::install_python) using this system

```
sessionInfo()
R version 4.4.0 (2024-04-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default


other attached packages:
 [1] reticulate_1.38.0        tensorflow_2.16.0.9000   keras_2.15.0             reshape2_1.4.4           dplyr_1.1.4             
 [6] gridGraphics_0.5-1       gridExtra_2.3            cowplot_1.1.3            xtable_1.8-4             synthpop_1.8-0          
[11] surrosurv_1.1.26         coxme_2.2-20             bdsmatrix_1.3-7          survival_3.5-8           metafor_4.6-0           
[16] numDeriv_2016.8-1.1      Matrix_1.7-0             meta_7.0-0               metadat_1.2-0            st4sas_1.0              
[21] maxLik_1.5-2.1           miscTools_0.6-28         gcipdr_0.0               ggplot2_3.5.1            mvtnorm_1.2-5           
[26] cubature_2.1.1           moments_0.14.1           JohnsonDistribution_0.24

```