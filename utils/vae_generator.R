# install the development version of packages, in case the
# issue is already fixed but not on CRAN yet.


if (Sys.info()['sysname'] == "Windows") {
  if (!require("keras") & !require("tensorflow") & !require("reticulate"))
  {
    install.packages("remotes")
    remotes::install_github(sprintf("rstudio/%s", c("reticulate", "tensorflow", "keras")))
    
    reticulate::install_python()
    keras::install_keras()
    
    # issues installing keras in windows. See : https://stackoverflow.com/questions/64774459/error-python-module-tensorflow-was-not-found-rstudio-windows10-path-problem
    
  }

} else {
  
  if (!require("reticulate")) {

    remotes::install_github("reticulate")
    reticulate::install_python()
    reticulate::virtualenv_create("r-venv", version = "3.10.14")
    
  }

}

if (Sys.info()['sysname'] != "Windows"){
  
  reticulate::use_virtualenv("~/.virtualenvs/r-keras/", required = TRUE)
  reticulate::py_config()
  
} else {
  
  reticulate::use_virtualenv("r-tensorflow")
  reticulate::py_config()
}

tensorflow::as_tensor("Hello World")

tensorflow::tf$constant("Hello TensorFlow!")

if (!require("keras")) {
  ## devtools::install_github("rstudio/keras")
  install.packages("keras")
  library(keras)
  ## install_keras()
}


if (!require("tensorflow")) {

  install.packages("tensorflow")
  library(tensorflow)
 # tensorflow::install_tensorflow(envname = "r-tensorflow")

}


# update: code appears to be broken with later package release. See https://divingintogeneticsandgenomics.com/post/how-to-code-a-variational-autoencoder-vae-in-r-using-mnist-dataset/ for running the network

### code based on: https://github.com/rstudio/keras/blob/master/vignettes/examples/variational_autoencoder.R

## see: https://cran.r-project.org/web/packages/keras/vignettes/faq.html

# With TF-2, you can still run this code due to the following line:
if (tensorflow::tf$executing_eagerly())
  tensorflow::tf$compat$v1$disable_eager_execution()


## use_session_with_seed(42)  ## this is broken under TF2, waiting to be fixed (https://github.com/rstudio/keras/issues/890)

tensorflow::tf$random$set_seed(0)  # temporary solution to make VAE results reproducible ...

K <- keras::backend()


# Parameters --------------------------------------------------------------

original_dim <- 5L
latent_dim <- 2L
intermediate_dim <- 10L
epsilon_std <- 1.0

# Model definition --------------------------------------------------------

x <- layer_input(shape = c(original_dim))
h <- layer_dense(x, intermediate_dim, activation = "relu")
z_mean <- layer_dense(h, latent_dim)
z_log_var <- layer_dense(h, latent_dim)

### reparametriztion trick --> z = mu + sigma*epsilon, epsilon~N(0,1)

sampling <- function(arg){
  z_mean <- arg[, 1:(latent_dim)]
  z_log_var <- arg[, (latent_dim + 1):(2 * latent_dim)]

  epsilon <- k_random_normal(
    shape = c(k_shape(z_mean)[[1]]),
    mean=0.,
    stddev=epsilon_std
  )

  z_mean + k_exp(z_log_var/2)*epsilon
}

# note that "output_shape" isn't necessary with the TensorFlow backend
z <- layer_concatenate(list(z_mean, z_log_var)) %>%
  layer_lambda(sampling)

# we instantiate these layers separately so as to reuse them later
decoder_h <- layer_dense(units = intermediate_dim, activation = "relu")
decoder_mean <- layer_dense(units = original_dim, activation = "relu")

## define decoding layers .... for the loss
h_decoded <- decoder_h(z)
x_decoded_mean <- decoder_mean(h_decoded)


## define decoding layers .... for the generator (note difference with above)
decoder_input <- layer_input(shape = latent_dim)
h_decoded_2 <- decoder_h(decoder_input)
x_decoded_mean_2 <- decoder_mean(h_decoded_2)


#### FULL NETWORK

# end-to-end autoencoder (e.g. from input to output)
vae <- keras_model(x, x_decoded_mean)

# encoder, from inputs to latent space
encoder <- keras_model(x, z_mean)

# generator, from latent space to reconstructed inputs
generator <- keras_model(decoder_input, x_decoded_mean_2)


vae_loss <- function(x, x_decoded_mean){
  xent_loss <- (original_dim/1.0)*loss_mean_squared_error(x, x_decoded_mean)
  kl_loss <- -0.5*k_mean(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var),
                         axis = -1L)
  xent_loss + kl_loss
}

vae %>%
  compile(
    optimizer = "rmsprop",
    loss = vae_loss)


#####################################################################
## vae wrapper
# argument 'repeat_data' clone input data with noise so many times as indicated (default = 100). This is akin to reinforcement learning to account for small sample size and create pseudo observations to artificially inflate the number of examples (data augmentation).

VAE <- function(x_train, H, epochs = 30L, repeat_data = 100, batch_split = 0.3)
{
    names <- colnames(x_train)
    x_train <- as.matrix(x_train)
    orig_dim <- dim(x_train)  # original dimension
    isbin <- apply(x_train, 2, is.binary)  # func 'is.binary' from gcipdr package
    x_train_augm <-  do.call(
        "rbind",
        lapply(1:repeat_data, function(i){  ## data augmentation
            out <- as.data.frame(matrix( nrow= orig_dim[1], ncol = orig_dim[2] ))
            colnames(out) <- names
            out[, isbin] <-  apply(as.matrix(x_train[, isbin]), 2, convert_binary)  ## smooth binary always
            if(i < 2)
                out[, !isbin] <- x_train[, !isbin]
            else   # perturb training data
                out[, !isbin] <- apply(as.matrix(x_train[, !isbin]), 2,
                                       function(x) x + rnorm(orig_dim[1], sd = 0.01))
            return(out)
        }
        )
    )
    examples_dim <- dim(x_train_augm)[1]   # dimension of augmented examples
    batch_size <- ifelse(examples_dim < 800, round(examples_dim*batch_split, 0), 100L)
    x_train_augm <- as.matrix(x_train_augm)


                                        # Model training ----------------------------------------------------------
    print(system.time(
        vae %>% fit(
                    x_train_augm, x_train_augm, ## x is the label of itself (autoencoding)
                    shuffle = TRUE,
                    epochs = epochs,
                    batch_size = batch_size,
                    validation_split = 0.2 # 20% of data used as validation for the loss
                )
    )
    )
### generate artificial data
    return(  lapply(1:H, function(i)
    {
        out <- predict(generator, matrix( MASS::mvrnorm(orig_dim[1], c(0,0), diag(2) ), ncol = 2) )
        out[, isbin] <- apply(out[, isbin], 2, function(x) ifelse(x < 0.5, 0, 1))  # reconvert to binary
        colnames(out) <- names
        return(out)
    }
    )
    )
}




############ DATA PREPROCESSING ###########
#### map binary into real line with knot point on zero


logistic <- function(x) 1 / (1 + exp(-x))

                                        # Two modifications in the VAE: a) convert binary digits, b) augment data to inflate number of examples.
                                        # a) we assign a fixed random Uniform quintile to 0 and 1. This seems to work better than creating a smooth line between 0 and 1, or between -5 and 5, where the latter must later be mapped onto [0,1] via sigmoid activation. Also it works better than just keeping original binary values. Data augmentation and a reasonable number of optimization epochs are needed for acceptable learning.


convert_binary <- function(x){

   if (!is.binary(x))
       stop("x must be binary")
    y <- x
    y <- ifelse( y < 1, runif(1, 0, 0.5), runif(1, 0.525, 1) )
    return(y)
}




