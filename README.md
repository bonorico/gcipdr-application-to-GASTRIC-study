# gcipdr-application-to-GASTRIC-study

Run file "main.R" for main analysis.

Run files "sim_study*.R" for simulation analyses (beware each file can take up to 24 hours to execute on 2 cores). Simulations were made with R 3.6.3 under a Linux 5.4.0-81-generic 18.04.1-Ubuntu system on a Intel(R) Core(TM) i7-2640M CPU @ 2.80GHz machine.  

This analysis was presented in [Survival analysis without sharing of individual patient data by using a Gaussian copula](https://onlinelibrary.wiley.com/doi/10.1002/pst.2415)

# Detail

Tensorflow installation for VAE might be problematic. If you are using Linux, in [vae_generator.R](https://github.com/bonorico/gcipdr-application-to-GASTRIC-study/blob/main/utils/vae_generator.R) be sure to have a fresh python installation and virtual environment:
```
reticulate::install_python()
reticulate::virtualenv_create("r-venv", version = "3.10.14")
```
If you are using Windows, you might want to follow [this](https://github.com/rstudio/tensorflow/issues/510#issuecomment-1045120710).

