## appendix_res.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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


## IPD KM curves

par(mfrow = c(1, 2))
 plot(survfit(fos, gastadj), fun = "cumhaz", lty = c(2, 1), main = "NA estimate")  # NA estimate  - log Surv
  plot(survfit(fos, gastadj), lty = c(2, 1), main = "KM estimate")

##### survival tables

kmOSmed <- summary(survfit(fos, data = gastadj))
    kmDFSmed <- summary(survfit(fdfs, data = gastadj))


OS1L <- kmOSmed$lower[kmOSmed$strata == "trt=1"] 
OS1Lt <- kmOSmed$time[kmOSmed$strata == "trt=1"]
OS1L <- OS1L[OS1Lt > 65 & OS1Lt <= 100]
OS1Lt <- OS1Lt[OS1Lt > 65 & OS1Lt <= 100]
OS0U <- kmOSmed$upper[kmOSmed$strata == "trt=0"] 
OS0Ut <- kmOSmed$time[kmOSmed$strata == "trt=0"]
OS0U <- OS0U[OS0Ut > 65 & OS0Ut <= 100]
OS0Ut <- OS0Ut[OS0Ut > 65 & OS0Ut <= 100]


DFS1L <- kmDFSmed$lower[kmDFSmed$strata == "trt=1"] 
DFS1Lt <- kmDFSmed$time[kmDFSmed$strata == "trt=1"]
DFS0U <- kmDFSmed$upper[kmDFSmed$strata == "trt=0"] 
DFS0Ut <- kmDFSmed$time[kmDFSmed$strata == "trt=0"]

DFS1La <- DFS1L[DFS1Lt > 5 & DFS1Lt <= 17]
DFS1Lta <- DFS1Lt[DFS1Lt > 5 & DFS1Lt <= 17]
DFS0Ua <- DFS0U[DFS0Ut > 5 & DFS0Ut <= 17]
DFS0Uta <- DFS0Ut[DFS0Ut > 5 & DFS0Ut <= 17]

DFS1L <- DFS1L[DFS1Lt > 55 & DFS1Lt <= 100]
DFS1Lt <- DFS1Lt[DFS1Lt > 55 & DFS1Lt <= 100]
DFS0U <- DFS0U[DFS0Ut > 55 & DFS0Ut <= 100]
DFS0Ut <- DFS0Ut[DFS0Ut > 55 & DFS0Ut <= 100]



## SUPPLEMENT FIGURE. ipd KM (mention: a cloglog plot would confirm PH assumption maybe not entirely met -- less effect at the beginning)
png(file.path(outfolder, "IPD_KMos.png"), width = 600, height = 1000)
par(mfrow = c(2, 1))
#   plot(survfit(fos, gastadj), conf.int = TRUE, lty = c(2, 1), ylab = "Kaplan-Meier estimate", xlab = "Time (months)")
plot(survfit(fos, gastadj), lty = c(2, 1), lwd = c(2, 2), ylab = "Kaplan-Meier estimate", xlab = "Time (months)", main = "OS", ylim = c(0.2, 1), cex.axis = 1.5, cex.lab = 1.5, cex = 1.5, cex.main = 1.5)
points(OS1Lt, OS1L, type = "S", col = "red", lwd = 2)
points(OS0Ut, OS0U, type = "S", col = "green", lty = 2, lwd = 2)
#   plot(survfit(fdfs, gastadj), conf.int = TRUE, lty = c(2, 1), ylab = "Kaplan-Meier estimate", xlab = "Time (months)")

plot(survfit(fdfs, gastadj), lty = c(2, 1), lwd = c(2, 2), ylab = "Kaplan-Meier estimate", xlab = "Time (months)", main = "DFS", ylim = c(0.2, 1), cex.axis = 1.5, cex.lab = 1.5, cex = 1.5, cex.main = 1.5)
points(DFS1Lt, DFS1L, type = "S", col = "red", lwd = 2)
points(DFS0Ut, DFS0U, type = "S", col = "green", lty = 2, lwd = 2)
points(DFS1Lta, DFS1La, type = "S", col = "red", lwd = 2)
points(DFS0Uta, DFS0Ua, type = "S", col = "green", lty = 2, lwd = 2)
legend("topright", c("Control", "Treatment", "Upper CI control", "Lower CI treatment"), lty = c(2, 1, 2, 1), lwd = c(2, 2, 2, 2), col = c("black", "black", "green", "red"), bty = "n")
dev.off()



  kmOS <- stat(gastadj, fos, KM  = T)
    kmDFS <- stat(gastadj, fdfs, KM  = T)



##### GRAPHICAL GOODNESS (pseudodata)

vertical.pooled <- do.call(
    "rbind",
    lapply(c("gamma","johnson"),
           function(j)
               data.frame(
                   do.call(
                       "rbind",
                       lapply(poolartd[[j]],
                              function(x)
                                  x
                              )
                   ),
                   DATA = j
               )
           )
)

vertical.gast <- rbind(vertical.pooled,
                       na.omit(
                           data.frame(
                               Return.IPD.design.matrix(gastadj[, -1],
                                                        fill = T
                                                        )[, -1],
                               trialID = gastadj$trialID,
                               DATA = "Orig. IPD"
                           )
                       )
                       ) 
vertical.gast$DATA <- as.factor(vertical.gast$DATA)
levels(vertical.gast$DATA)[1:2] <- c("GC few. mom.", "GC all mom.")
vertical.gast$trt <- ifelse(vertical.gast$trt == 1, "Treatment", "Control" )
vertical.gast$statusT <- ifelse(vertical.gast$statusT == 1, "Deceased", "Alive" )
vertical.gast$statusS <- ifelse(vertical.gast$statusS == 1, "Relapsed", "Non-relapsed" )

#### histograms

plotdat <- vertical.gast[, c("timeS", "timeT","statusS", "statusT", "trt", "DATA")]

detail <- textsize <- theme(axis.text = element_text(size=18), axis.title = element_text(size=20), strip.text.x = element_text(size = 19), strip.text.y = element_text(size = 19), legend.title=element_blank(), legend.text=element_text(size=18), legend.position="bottom" )

mpT <- aes( x= timeT, y=..density.., fill = DATA, color = DATA)
plT <- ggplot(plotdat, mpT )
histT <- geom_density(alpha=0.2, position="identity", show.legend = FALSE)
panelT <- facet_grid( as.factor(trt) ~ as.factor(statusT) )
yl <- ylab("Empirical density")
xlT <- xlab("Time to death (months)")

 hist_timeT <- plT + histT + panelT + xlT + yl + theme_bw() + xlim(0, 150) + detail

#

mpS <- aes( x= timeS, y=..density.., fill = DATA, color = DATA)
plS <- ggplot( plotdat, mpS )
panelS <- facet_grid(as.factor(trt) ~ as.factor(statusS) )
histS <- geom_density(alpha=0.2, position="identity")
 xlS <- xlab("Time to relapse (months)")

 hist_timeS <- plS + histS + panelS + xlS + yl + theme_bw() + theme(legend.position="bottom" ) + xlim(0, 150) + detail


### 

png(file.path(outfolder, "hist_gast.png"), width = 800, height = 1000)

 plot_grid( hist_timeT, hist_timeS, ncol = 1, labels = "AUTO" )

dev.off()


### JOINT DEPENDENCE

aesS <- aes(x = timeS, y = timeT, fill = DATA, color = DATA)
aesT <- aes(x = timeT, y = timeS, fill = DATA, color = DATA)

plS <- ggplot(plotdat, aesS)
plT <- ggplot(plotdat, aesT)

pointsT <- geom_point( alpha = 0.2, show.legend = FALSE )
pointsS <- geom_point( alpha = 0.2) 

jT <- plT + pointsT + panelT + theme_bw() + xlS + ylab("Time to death (months)") + detail + xlim(0, 600) + ylim(0, 600)

jS <- plS + pointsS + panelS + theme_bw() + ylab("Time to relapse (months)") + xlT + guides(colour = guide_legend(override.aes = list(alpha = 0.8, size = 4))) + detail + xlim(0, 600) + ylim(0, 600)

#
png(file.path(outfolder, "joint_gast.png"), width = 800, height = 800)

 plot_grid( jT, jS, ncol = 1, labels = "AUTO" )

dev.off()





