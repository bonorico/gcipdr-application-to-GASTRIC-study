## my_utils.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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


ztest <- function(d) rbind(d[1]/d[2], pnorm(d[1]/d[2], lower = TRUE)*2)

Ztest <- function(d) apply(d, 2, function(x) ztest(x))

ZtestS <- function(d, drop = c(1,2))
{    
    n <- dim(d)[1]
    s <- which(d$Est. == "sd")
    inp <- as.matrix(d[c(s-1, s), -drop])
    out <- data.frame(Out. = unique(d[, 1]),
                      Est. = c("z.1", "p.1"),
                      Ztest(inp)
                      )
    res <- rbind(d[1:s, ], out, d[(s+1):n, ])    
    return(res)
}


ZtestL <- function(d,
                   keep = c("Death", "Disease"),
                   drop = c(1:2)
                   )
    do.call(
        "rbind",
        lapply(keep,
               function(o)
                   ZtestS(d[d$Out. == o, ],
                          drop
                          )
               )
    )


#


LRtest <- function(LRnew, LRref, DFnew, DFref)
{
    lr <- LRnew-LRref ## modelref must be smaller than model new
    df <- DFnew-DFref
    p <- pchisq( lr, df, lower = FALSE) 

    data.frame(Chi = lr, df = df, p = p)
} 




                                        # small functions to boldface results if condition is met

Boldnum <- function(x)
{
    new <- ifelse(is.na(x) | is.null(x),
                  NA,
                  paste0("\\textbf{", x, "}")
                  )
    return(new)
}

                                        #

Boldnums <- function(y, tresh,
                     err, roundup = 2)
{
                                        # this function evaluates err >= tresh as condition
    y <- as.matrix(y)
    out <- ifelse(err >= tresh,
                  Boldnum(round(y, roundup) ),
                  round(y, roundup)
                  ) 
    
    return(out)

}

                                        #


                                        #
rsq <- function(y, fit, p = 1, adj = FALSE)
{
    x <- na.omit(cbind(y, fit))
    n <- dim(x)[1]
    k <- 1
    sst <- var(as.numeric(x[, 1]))*(n-1)              # total sum of squares
    ssr <- sum((as.numeric(x[, 1]) - as.numeric(x[, 2]))^2)   # residual sum of squares
    if (adj)
        k <- (n-1)/(n-p-1)
        
    rsq <- 1-(ssr/sst)*k
    mse <- ssr/n
    return(c(Rsq = rsq, RMSE = sqrt(mse)))

}

                                        #
print_synopsis <- function(X, names_X = NULL){
    # print a relationship synopsis (pearson cor, R-sq, RMSE) of P variables in X (all pairs).
    # Args: X -> numerical matrix
              
    P <- dim(X)[2]
    K <- P*(P-1)/2
    
    if (is.null(names(X)) & is.null(names_X))
        names(X) <- paste0("V", 1:P)
    if (!is.null(names_X))
        names(X) <- names_X

    out <- do.call("rbind",
                   lapply(1:K, function(k)
                   {
                       combos <- combn(P, 2)
                       i <- combos[1, k]
                       j <- combos[2, k]

                       cor <- cor.test(X[, i],
                                       X[, j],
                                       method = "pearson")  # TODO: can be set as Arg.
                       oth <- rsq(X[, i],
                                  X[, j])
                       Rsq <- oth[1]
                       rmse <- oth[2]
                       
                       data.frame(
                           pair = paste0(
                               "(",names(X)[i], ") vs (",
                               names(X)[j],")"),
                           cor = cor$estimate,
                           corpval = cor$p.value,
                           Rsq = Rsq,
                           RMSE = rmse

                       )
                       
                   }  ) )
    return(out)

}

#

km_ks_test <- function(Xlong, t_value = "time", prob_value = "km", from = quote(Estimate))
{
                                        # Xlong = matrix, t_value -> numeirc, prob_value -> numeirc, from = factor
    groups <- levels(Xlong[[from]])

    P <- length(groups)
    K <- P*(P-1)/2
#    browser()
    out <- do.call("rbind",
                   lapply(1:K, function(k)
                   {
                       combos <- combn(P, 2)
                       i <- combos[1, k]
                       j <- combos[2, k]
                       
                       X <- subset(Xlong, eval(from) == groups[i], select = c(t_value, prob_value))
                       Y <- subset(Xlong, eval(from) == groups[j], select = c(t_value, prob_value))
                       
                       x <- X[, t_value]
                       y <- Y[, t_value]
                       Fx <- X[, prob_value]
                       Fy <- Y[, prob_value]
                       
                       Sdiff <- ks_test(x, y, Fx, Fy)    # Kolmogorov-Smirnov test (unweighted !!) --> just a crude difference
                       
                       data.frame(
                           pair = paste0(
                               "(",groups[i], ") vs (",
                               groups[j],")"),
                           KSstat = Sdiff["D"],
                           level = Sdiff["level"],
                           refused = Sdiff["H1"]

                       )
                       
                   }  ) )
    return(out)
    
}

                                        #
filter_time <- function(x, y)
{
    # y = reference
    # returns rank(x): x-y = min
    ty <- which(y >= min(x))
    y <- y[ty] 
    if (length(y) > length(x))
        warning("len x must >= len y")
    rx <- range(x)
    ry <- range(y)
    u <- c(rx[2], ry[2])
    l <- c(rx[1], ry[1])
    ovrc <- (rx[1] < ry[2]) | (ry[1] < rx[2])   # overlap condition
    dx <- diff(rx)
    dy <- diff(ry)
    ovr <- min(c(dx, dy))/(max(dx, dy))  # ranges overlap amount, provided they do overlap
    if (ovrc & (ovr < 0.5))
        message("x and y ranges have less than 50% overlap")    
    if (!ovrc)
        warning("x and y ranges do not overlap")
      
    tx <- unlist(
        lapply(y,
               function(i)
               {
                   j <- which(x <= i)
                   jmin <- tail(x[j], n = 1)
                   rank <- tail(which(x == jmin), n = 1)                   
                   return(rank)
               }
                                             
               )
    )
    return(list(tx, ty))
}
                                        #

ks_test <- function(x, y, Fx, Fy, alpha = 0.05, conservative = FALSE)
{
    x <- sort(x)
    y <- sort(y)
    Fx <- Fx[order(x)]
    Fy <- Fy[order(y)]
    
    nx <- length(x)
    ny <- length(y)
    if (nx >= ny)
    {
        cat("Using y as reference\n")
        xy <- filter_time(x, y)
        Fx <- Fx[xy[[1]]]
        Fy <- Fy[xy[[2]]]
    }
    else
    {
        cat("Using x as reference\n")
        yx <- filter_time(y, x)
        Fy <- Fy[yx[[1]]]
        Fx <- Fx[yx[[2]]]
    }

    n <- length(Fx)
    m <- length(Fy)
    if (n != m)
        stop ("Fx Fy differ in length")

    ca <- -log(alpha/2)*(1/2)
    w <- (n + m)/(n*m)    # equals 2/n --> based on actually used points
    if (conservative)
        w <- (nx + ny)/(nx*ny)  # based on original (maybe unused) sample size
    level <- sqrt(ca*w)
    diffs <- abs(Fx - Fy)
    D <- max(diffs)[1]    # KS statistic = supremum_t over abs differences (drop ties)
    H1 <- (D > level)    # test -> refuse H0 if TRUE
    if (H1)
        cat("H0 rejected: distributions are different\n")
    else
        cat("H0 accepted: distributions might be similar\n")
        
    return(c(D = D, level = level, H1 = H1))
                        
}




############ other plottin tools ... contours


if (!require("reshape2")) {
    install.packages("reshape2")
    library(reshape2)
}


contour.data <- function(D, cname1, cname2,
                         breaks = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99),
                         cumulative = TRUE)
{
    kk <- MASS::kde2d(D[, cname1],
                      D[, cname2],
                      n = 100
                      )  

    dx <- diff(kk$x[1:2])
    dy <- diff(kk$y[1:2])

    sz <- sort(kk$z)
    c1 <- cumsum(sz) * dx * dy
    
    dimnames(kk$z) <- list(kk$x,kk$y)
    dc <- melt(kk$z)    
    colnames(dc)[1:2] <- c(cname1, cname2)

    dcs <- dc[order(dc$value), ]  ## order grid along cumulative sum

    Cdat <- data.frame(dcs, cumul.dens = 1-c1)

    contourbreaks <- NULL

    if (!cumulative)
        contourbreaks <- do.call(
            "rbind",
            lapply(breaks,
                   function(b)
                       Cdat[
                           Cdat$cumul.dens == max(
                                                  Cdat$cumul.dens[Cdat$cumul.dens <= b]
                                              ), ]
                   )
        )

    return(list(contours= Cdat,
                breaks = contourbreaks
                )
           )

}
                           



#### summarizing tools


summ <- function(x,
                 mean = FALSE,
                 onlymean = FALSE,
                 label = FALSE,
                 rd = 1)
    ifelse(label,
    ifelse(mean & !onlymean,
           paste0(round(mean(x, na.rm = T), rd),
                  " (SD = ", round(sd(x, na.rm = T), rd),
                  "; range = ",
                  round(min(x, na.rm = T), rd),
                  ", ",
                  round(max(x, na.rm = T), rd), ")" ),
    ifelse(mean & onlymean,
           round(mean(x, na.rm = T), rd),
           paste0("(SD = ", round(sd(x, na.rm = T), rd),
                  "; range = ", round(min(x, na.rm = T), rd),
                  ", ",
                  round(max(x, na.rm = T), rd), ")" ) )
    ),
    ifelse(mean & !onlymean,
           paste0(round(mean(x, na.rm = T), rd),
                  " (", round(sd(x, na.rm = T), rd),
                  "; ", round(min(x, na.rm = T), rd),
                  "/", round(max(x, na.rm = T), rd),
                  ")" ),
    ifelse(mean & onlymean,
           round(mean(x, na.rm = T), rd),
           paste0("(",
                  round(sd(x, na.rm = T), rd),
                  "; ", round(min(x, na.rm = T), rd),
                  "/", round(max(x, na.rm = T), rd), ")" ) )
    )
    )
        
#

corr_names <- function(names, lowertri = TRUE)
{
    out <- matrix(
        unlist(
            lapply(names,
                   function(i)
                       lapply(names,
                              function(j)
                                  paste(i,j, sep = "-")
                              )
                   )
        ),
        nrow = length(names)
    )
    if (lowertri)
        out <- out[lower.tri(out)]
    return(out)
    
} 

#

summ_summs <- function(sumlist,
                       mean = FALSE,
                       onlymean = FALSE,
                       label = FALSE, rd = 1)
{

    varnames <- rownames(sumlist[[1]]$moms)
    momnames <- colnames(sumlist[[1]]$moms)
    N <- length(sumlist)
    p <- length(varnames); P <- p*(p-1)/2

    momsum <-  do.call(
        "cbind",
        lapply(momnames,
               function(m)
                   do.call(
                       "rbind",
                       lapply(varnames,
                              function(n)
                                  data.frame(out = summ(
                                                 unlist(
                                                     lapply(1:N,
                                                            function(i)
                                                                sumlist[[i]]$moms[n, m]
                                                            )
                                                 ),
                                                 mean,
                                                 onlymean,
                                                 label,
                                                 rd
                                             )
                                             )
                              )
                   )
               )
    )
    colnames(momsum) <- momnames

    corsum <- do.call(
        "rbind",
        lapply(1:P,
               function(j)
                   data.frame(cor = summ(
                                  unlist(
                                      lapply(1:N,
                                             function(i)
                                                 sumlist[[i]]$corr[lower.tri(sumlist[[i]]$corr)][j]
                                             )
                                  ),
                                  mean,
                                  onlymean,
                                  label,
                                  rd
                              )
                              )
               )
    )

    res <- rbind(
        data.frame(momsum,
                   cor = NA),
        data.frame(mx = NA,
                   sdx = NA,
                   skx = NA,
                   ktx = NA ,
                   corsum
                   )
    )

    rownames(res) <- c(varnames, corr_names(varnames))
    
    return(res)

}  



############## diagnostics on artif data


diagnostic_summ <- function(simdlong, drop)
{
    istsimulsumm <- lapply(simdlong,
                           function(j)
                           {
                               data <- j[, -drop]
                               corr <- momcor(apply(data, 2, as.numeric))
                               list(mx= apply(data, 2, mean, na.rm = T),
                                    sdx = apply(data, 2, sd, na.rm = T),
                                    skx = apply(data, 2, skewness, na.rm = T),
                                    ktx = apply(data, 2, kurtosis, na.rm = T),
                                    corr = corr[lower.tri(corr)]
                                    )
                           }
                           )
    aversumm <- lapply(c("mx", "sdx", "skx", "ktx", "corr"),
                       function(s)
                           apply(
                               do.call(
                                   "rbind",
                                   lapply(istsimulsumm ,
                                          function(j)
                                              j[[s]]
                                          )
                               ), 2, mean, na.rm = T
                           )
                       )
    names(aversumm[[5]]) <- corr_names(names(aversumm[[1]]))

    names(aversumm) <- c("mx", "sdx", "skx", "ktx", "corr")

    return(aversumm)

}



## empirical computation of time-averaged HR 

time_aver_HR <- function(ipd, b1, b2)
{
                                        #      ipd <- ipd[ipd$status == 1, ]  # to include/exclude censoring is irrelevant due to independence
    t1 <-  max(ipd$time[ipd$scen == 1], na.rm = T) # Weibull time 1
    t2 <-  max(ipd$time[ipd$scen == 2], na.rm = T) # Weibull time 2
    return(  (b1*t1 + b2*t2)/(t1+t2)   )
}





### # NOT needed
pivot_wide <- function(d, target_value = "km", from = "Estimate", length_from = "Orig. IPD")
{
    
    ## sel <- unlist(
    ##     lapply(c(target, other),
    ##            function(j)
    ##                which(names(d) == j))
    ## )
    sel <- which(names(d) == target_value)
    ref <- d[[from]]
    index <- levels(ref)
    grab <- lapply(index,
                   function(j)
                       d[which(ref == j), sel]
                   )
    
    N <- length(
        grab[[which(index == length_from)]]
        )

    out <- do.call("cbind",
                   lapply(grab,
                          function(j)
                          {
                              n <- length(j)
                              if(n >= N)
                                  return(j[1:N])
                              else
                                  return(
                                      c(j,
                                        rep(NA, (N-n))
                                        )
                                  )

                          }


                          )
                   )

    
    colnames(out) <- index
    return(out)           
}
