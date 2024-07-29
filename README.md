# gcipdr-application-to-GASTRIC-study

Run file "main.R" for main analysis.

Run files "sim_study*.R" for simulation analyses (beware each file can take up to 24 hours to execute on 2 cores). Simulations were made with R 3.6.3 under a Linux 5.4.0-81-generic 18.04.1-Ubuntu system on a Intel(R) Core(TM) i7-2640M CPU @ 2.80GHz machine.  

This analysis was presented in [Survival analysis without sharing of individual patient data by using a Gaussian copula](https://onlinelibrary.wiley.com/doi/10.1002/pst.2415)

# Detail

Tensorflow installation for VAE might be problematic. The code attempts to distinguish between OS (Windows or not) and tailor installation accordingly. Because of rapid development of packages tensorflow and keras the current VAE code might become outdated any time soon. Efforts to mantain the VAE code cannot be guaranteed.

