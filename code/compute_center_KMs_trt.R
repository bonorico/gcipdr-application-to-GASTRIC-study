## compute_center_KMs.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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


##### CENTER-EFFECT ON SURVIVAL UNDER TREATMENT


#### TWO STAGE (parametric) APPROACH

trthazdata <- subset(gastadj, trt == 1)

kmtrtma <- do.call("rbind", lapply( c("Death", "Disease"), function(i) do.call("rbind", lapply(trialnames, function(j){
    out <- with(subset(trthazdata, trialID == j), subset(kmparam(
                                                       switch(i, "Death" = timeT, "Disease" = timeS),
                                                       switch(i, "Death" = statusT, "Disease" = statusS), trt, 300, NULL), trt == 1))
        data.frame(Outcome = i, Estimate = "Meta an.", t = out$time, km = out$km, Center = j)
     } ) ) ) )


#### ONE STAGE with PSEUDODATA
mod_km <- c("Death" = Surv(timeT, statusT) ~ 1,
            "Disease" = Surv(timeS, statusS) ~ 1)

system.time(
    RES_model_kmtrt <- lapply(mod_km, function(m)
        lapply(pd, function(d) mclapply(1:length(trialnames), function(j)
            BootStat(
                lapply(d, function(x) subset(x, trt == 1 & trialID == as.character(j) )),
                "km", reg_formula = m) ## km

            )   
            )    )
)

#
kmtrtdata <- do.call("rbind", lapply(RES_model_kmtrt, function(m)
    do.call("rbind", lapply(m, function(d) do.call("rbind", lapply(1:length(trialnames), function(j)
        data.frame(t = d[[j]]$Bag$time, km = d[[j]]$Bag$km, Center = as.character(j))        
                                                                   ) ) )) ) )

kmtrtdata <- data.frame(
    Outcome = as.factor(
        unlist(lapply(rownames(kmtrtdata), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1] ) ) ),
    Estimate = as.factor(unlist(lapply(rownames(kmtrtdata), function(x) strsplit(x, ".", fixed = TRUE)[[1]][2] ) ) ),
    kmtrtdata
        )

                                        #
kmtrtdata <- rbind( kmtrtdata, kmtrtma)

levels(kmtrtdata$Estimate)[2:6] <- c("Dist. boot.", "GC all mom.", "GC few. mom.", "Orig. IPD", "IPD boot.")
kmtrtdata$Estimate <- factor(kmtrtdata$Estimate, levels = c("Orig. IPD", "IPD boot.", "Dist. boot.", "CART",
                                                      "VAE", "GC few. mom.", "GC all mom.", "Meta an."))

ordc <- as.character(sort(as.numeric(kmtrtdata$Center)))
kmtrtdata$Center <- factor(kmtrtdata$Center, levels = unique(ordc))

 ## plot center-specific log Baseline hazards (log-log baseline KMs)

kmtrt1 <- subset(kmtrtdata, Outcome == "Death" & (Estimate == "GC all mom." | Estimate == "Orig. IPD"))


km_annot <- function(dat, centers = quote(Center == "2" | Center == "3" | Center == "4" | Center == "8"))
{
    dat %>%
        group_by(Center, Estimate) %>%
        summarise(t = last(t) + 4, km = last(km)) 
}


center_label <-  subset(km_annot(kmtrt1), Center == "2" | Center == "3" | Center == "4" | Center == "8" )

    
mp <- aes( x= t, y= km, color = Center)
kmtrtplot <- ggplot( kmtrt1 , mp )   # plot only OS (no new insight in DFS)
steps <- geom_step(alpha=0.5, size = 1.3)
panel <- facet_grid( ~ as.factor(Estimate) )
xl <- xlab("Time (months)");  yl <- ylab("Kaplan-Meier")
textsize <- theme(axis.text = element_text(size=18), axis.title = element_text(size=20), strip.text.x = element_text(size = 19), legend.title=element_text(size=18), legend.text=element_text(size=16), legend.position="bottom" )    
annot <- geom_text(aes(label = Center), data = center_label,
                   show.legend = FALSE,
                   colour = "black",
                   size = 6 )


## pdf(file.path(outfolder, "KM_trt_centers_plot_death_John.eps"), width = 12.4, height = 8)

main3 <- kmtrtplot + steps + panel + xl + yl + theme_bw() + textsize + annot

##dev.off()

panel <- facet_wrap( ~ as.factor(Estimate), nrow = 2 )
kmtrtplot_dea <- ggplot( subset(kmtrtdata, Outcome == "Death") , mp )   # plot only OS (no new insight in DFS)
center_label_all <-  subset(km_annot(kmtrtdata),
                            Center == "2" | Center == "3" | Center == "4" | Center == "8" | Center == "7" | Center == "9" )
center_label_all[center_label_all$Estimate == "CART" & center_label_all$Center == "9", "km"] <- 0.35
center_label_all[center_label_all$Estimate == "CART" & center_label_all$Center == "9", "t"] <- 135

annot_all <- geom_text(aes(label = Center), data = center_label_all,
                   show.legend = FALSE,
                   colour = "black",
                   size = 6 )


pdf(file.path(outfolder, "KM_trt_centers_plot_death.pdf"), width = 25, height = 16)

 kmtrtplot_dea + steps + panel + xl + yl + theme_bw() + textsize + xlim(0, 350) + annot_all  

dev.off()


print("New Objects: 'main3'")

                                        #

#### NO NEW INSIGHT: DON'T PLOT
kmtrtplot_dis <- ggplot( subset(kmtrtdata , Outcome == "Disease") , mp )   

#### pdf("BH_centers_plot_disease.eps", width = 25, height = 8)

 kmtrtplot_dis + steps + panel + xl + yl + theme_bw() + textsize + xlim(0, 300)  + annot

#### dev.off()


### CLARIFY OUTLIER IN CENTER 2

out <- survfit(Surv(timeT, statusT)~trt, subset(gastadj, trialID == 2 ))
 plot(out, lty = c(2, 1), main = "KM estimate")

summary(out)
