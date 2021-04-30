## subgroup_analysis.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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

#### SUB-GROUP ANALYSIS
#### NON-LINEARITY: PH assumption violated ? 

subga <- function(simd, formula, tcutoff)
{
    within <- BootStat(
        lapply(simd,
               function(x)
                   subset(x,
                          subset = timeT <= tcutoff
                          )
               ),
        "fixed", formula
    )
                                        # HR within cutcoff
    over <- BootStat(
        lapply(simd,
               function(x)
                   subset(x,
                          subset = timeT > tcutoff
                          )
               ),
        "fixed", formula
    )
                                        # HR over cutoff

    out <-  data.frame(within$Bag$coef,
                       within$Bag$sd,
                       over$Bag$coef,
                       over$Bag$sd
                       )

    colnames(out) <- c(
        paste0("HR_t<=", tcutoff),
        paste0("HRsd_t<=", tcutoff),
        paste0("HR_t>", tcutoff),
        paste0("HRsd_t>", tcutoff)
    )

    return(out)
}




plot(
    survfit(fos, gastadj),
    fun = "cloglog",
    lty = c(2, 1),
    ylab = "Cloglog Survival",
    xlab = "Scaled time ",
    xlim = c(18, 400),
    ylim = c(-1.5, 0.6)
)
legend("bottomright",
       c("Control", "Treatment"),
       lty = c(2, 1),
       bty = "n"
       )

                                        #

table2_bold <- table2 <- do.call("rbind", lapply( pd, function(x) subga(x, fos_1a, 25) ) )

rownames(table2) <- c("Orig. IPD", "IPD boot.", "Dist. boot.", "CART", "VAE", "GC few. mom.", "GC all mom.")


table2_bold[, 1] <- Boldnums(
    table2[, 1],
    abs(qnorm(0.01/2)),
    abs(table2[, 1]/table2[, 2]),
    round_up
)

table2_bold[, 3] <- Boldnums(
    table2[, 3],
    abs(qnorm(0.01/2)),
    abs(table2[, 3]/table2[, 4]),
    round_up
)

table2_bold[, c(2,4)] <-  round(table2[, c(2,4)], round_up)

print("New Objects: 'table2', 'table2_bold'")

#### BIG LIMITATION (FUTURE WORK): the linear correlation assumption misses structure in the data. We need a concise statistic to consicely report data non-linearity. A concise expression y~x is a good candidate. Future work: Start from gastadj example statusT~timeT to propose such concise expression (a logistic piecewise report ?)
