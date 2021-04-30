## compute_center_BHs.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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




###### INSPECT ORIGINAL BASELINE HAZARDS in each center: Have they the same shape ? Which RE model is more appropriate ? 
### relates to MODEL_1d, e, g


#### TWO STAGE (parametric) APPROACH

basehazdata <- subset(gastadj, trt == 0)

bhma <- do.call("rbind", lapply( c("Death", "Disease"), function(i) do.call("rbind", lapply(trialnames, function(j){
    out <- with(subset(basehazdata, trialID == j), subset(kmparam(
                                                       switch(i, "Death" = timeT, "Disease" = timeS),
                                                       switch(i, "Death" = statusT, "Disease" = statusS), trt, 300, NULL), trt == 0))
        data.frame(Outcome = i, Estimate = "Meta an.", t = out$time, bh = -log(out$km), Center = j)
     } ) ) ) )


#### ONE STAGE with PSEUDODATA

model_bh <- c("Death" = Surv( timeT, statusT )~1,
              "Disease" = Surv( timeS, statusS )~1)
#### RUN MODEL

 system.time(
RES_model_bh <- lapply(model_bh, function(m)
    lapply(pd, function(d) mclapply(1:length(trialnames), function(j)
        BootStat(
            lapply(d, function(x) subset(x, trt == 0 & trialID == as.character(j) )),
            "km", reg_formula = m, cumhaz = TRUE) ## km-like --> nelson-aalen

                                  )   
           )    )
)

#
bhdata <- do.call("rbind", lapply(RES_model_bh, function(m)
    do.call("rbind", lapply(m, function(d) do.call("rbind", lapply(1:length(trialnames), function(j)
        data.frame(t = d[[j]]$Bag$time, bh = d[[j]]$Bag$km, Center = as.character(j))        
                                                                   ) ) )) ) )

bhdata <- data.frame(
    Outcome = as.factor(
        unlist(lapply(rownames(bhdata), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1] ) ) ),
    Estimate = as.factor(unlist(lapply(rownames(bhdata), function(x) strsplit(x, ".", fixed = TRUE)[[1]][2] ) ) ),
    bhdata
        )

                                        #
bhdata <- rbind( bhdata, bhma)

levels(bhdata$Estimate)[2:6] <- c("Dist. boot.", "GC all mom.", "GC few. mom.", "Orig. IPD", "IPD boot.")
bhdata$Estimate <- factor(bhdata$Estimate, levels = c("Orig. IPD", "IPD boot.", "Dist. boot.", "CART",
                                                      "VAE", "GC few. mom.", "GC all mom.", "Meta an."))

ordc <- as.character(sort(as.numeric(bhdata$Center)))
bhdata$Center <- factor(bhdata$Center, levels = unique(ordc))

bhdata$bh <- log( bhdata$bh )  # ... log hazard


bh_annot <- function(dat)
{
    dat %>%
    group_by(Center, Estimate) %>%
    summarise(t = last(t) + 4, bh = last(bh)) %>%
    subset(Center == "3" | Center == "4" | Center == "9" | Center == "13")
}

## plot center-specific log Baseline hazards (log-log baseline KMs)

bshaz1 <- subset(bhdata, Outcome == "Death" & (Estimate == "GC all mom." | Estimate == "Orig. IPD"))

center_label <-  bh_annot(bshaz1)

 mp <- aes( x= t, y= bh, color = Center)
bhkmplot <- ggplot( bshaz1 , mp )   # plot only OS (no new insight in DFS)
steps <- geom_step(alpha=0.5, size = 1.3)
panel <- facet_grid( ~ as.factor(Estimate) )
xl <- xlab("Time (months)");  yl <- ylab("Log Nelson-Aalen")
ylims <- ylim(-5, 1.5)
textsize <- theme(axis.text = element_text(size=18), axis.title = element_text(size=20), strip.text.x = element_text(size = 19), legend.position = "none")   ## legend.title=element_text(size=18), legend.text=element_text(size=16), 
annot <- geom_text(aes(label = Center), data = center_label,
                   show.legend = FALSE,
                   colour = "black", check_overlap = TRUE,
                   size = 6 )


##pdf(file.path(outfolder, "BH_centers_plot_death_John.eps"), width = 12.4, height = 8)

main4 <- bhkmplot + steps + panel + xl + yl + theme_bw() + textsize + ylims + annot
 
##dev.off()



panel <- facet_wrap( ~ as.factor(Estimate), nrow = 2 )
bhplot_dea <- ggplot( subset(bhdata, Outcome == "Death") , mp )   # plot only OS (no new insight in DFS)
center_label_all <-  bh_annot(bhdata)
annot_all <- geom_text(aes(label = Center), data = center_label_all,
                   show.legend = FALSE,
                   colour = "black",
                   size = 6 )


 pdf(file.path(outfolder, "BH_centers_plot_death.eps"), width = 25, height = 16)

 bhplot_dea + steps + panel + xl + yl + theme_bw() + textsize + ylims + xlim(0, 350) + annot_all

 dev.off()


print("New Objects: 'main4'")#

#### NO NEW INSIGHT: DON'T PLOT
bhplot_dis <- ggplot( subset(bhdata, Outcome == "Disease") , mp )   

#### pdf("BH_centers_plot_disease.eps", width = 25, height = 8)

bhplot_dis + steps + panel + xl + yl + theme_bw() + textsize + ylims + xlim(0, 300)  

#### dev.off()


