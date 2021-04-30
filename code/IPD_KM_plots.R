## IPD_KM_plots.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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


#### plot estimates from IPD_KM.R
## uses: kmdata


###  TREATMENT EFFICACY ON SURVIVAL  
## plot KM (overlapped to IPD)

## OS only

KM_GC_main <- subset(kmdata, Outcome == "Overall survival probability" & (Estimate == "Orig. IPD" | Estimate == "GC all mom.") )

mp0 <- aes( x= time, y= km, color = Treatment)
kmplot1 <- ggplot( KM_GC_main, mp0 )
steps <- geom_step(alpha=0.5, size = 1.3)
panel <- facet_wrap(~ as.factor(Estimate), ncol=2)
xl <- xlab("Time (months)");  yl <- ylab("Kaplan-Meier")
xlims <- xlim(0, 310)
textsize <- theme(axis.text = element_text(size=18), axis.title = element_text(size=20), strip.text.x = element_text(size = 19), legend.title=element_blank(), legend.text=element_text(size=18), legend.position="bottom" )


mp <- aes( x= time, y= km, color = Estimate)
panel2 <- facet_wrap(~ as.factor(Treatment), ncol=2)
kmplot2 <- ggplot( KM_GC_main, mp )

 main1 <-   kmplot1 + steps + panel + xl + yl + xlims + theme_bw() + textsize
 main2 <-   kmplot2 + steps + panel2 + xl + yl + xlims + theme_bw() + textsize


### FIGURE 1 OF MAIN MANUSCRIPT (OS Only)

pdf(file.path(outfolder, "km_plot_OS.eps"), width = 10.4, height = 12.4)

plot_grid(main1, main2, labels="AUTO", nrow = 2)

dev.off()


### KM plots overlapped (complete)

KM_MA <- subset(kmdata, Estimate == "Orig. IPD" | Estimate == "Meta an.")
KM_GC_gam <- subset(kmdata, Estimate == "Orig. IPD" | Estimate == "GC few. mom.")
KM_GC_john <- subset(kmdata, Outcome == "Disease free survival probability" &
                             ( Estimate == "Orig. IPD" | Estimate == "GC all mom."))
KM_CART <- subset(kmdata, Estimate == "Orig. IPD" | Estimate == "CART")
KM_BOOT <- subset(kmdata, Estimate == "Orig. IPD" | Estimate == "IPD boot.")
KM_DISTBOOT <- subset(kmdata, Estimate == "Orig. IPD" | Estimate == "Dist. boot.")
KM_VAE <- subset(kmdata, Estimate == "Orig. IPD" | Estimate == "VAE")

# 
panel <- facet_wrap(~ as.factor(Outcome) + as.factor(Treatment), ncol=2)

#
png(file.path(outfolder, "km_plot_GCgamma.png"), width = 1000, height = 800)

ggplot( KM_GC_gam , mp ) + steps + panel + xl + yl + xlims + theme_bw() + textsize

dev.off()

                                        #
png(file.path(outfolder, "km_plot_GCjohn.png"), width = 1000, height = 400)

ggplot( KM_GC_john , mp ) + steps + panel + xl + yl + xlims + theme_bw() + textsize

dev.off()


### FIG KM.synthpop supplement material

 png(file.path(outfolder, "km_plot_synpd.png"), width = 1000, height = 800)

ggplot( KM_CART, mp ) + steps + panel + xl + yl + xlims + theme_bw() + textsize 

dev.off()


### FIG Boot.KM supplement material

 png(file.path(outfolder, "km_plot_boot.png"), width = 1000, height = 800)

ggplot( KM_BOOT, mp ) + steps + panel + xl + yl + xlims + theme_bw() + textsize 

dev.off()


### FIG BootDist.KM supplement material

 png(file.path(outfolder, "km_plot_bootDist.png"), width = 1000, height = 800)

ggplot( KM_DISTBOOT , mp ) + steps + panel + xl + yl + xlims + theme_bw() + textsize 

dev.off()


### FIG KM.VAE supplement material

 png(file.path(outfolder, "km_plot_vae.png"), width = 1000, height = 800)

ggplot( KM_VAE , mp ) + steps + panel + xl + yl + xlims + theme_bw() + textsize 

dev.off()


### FIG KM.metaanal supplement material

 png(file.path(outfolder, "km_plot_meta.png"), width = 1000, height = 800)

ggplot( KM_MA , mp ) + steps + panel + xl + yl + xlims + theme_bw() + textsize 

dev.off()


#### KM ribbon plots only CIs (side-by-side to IPD)

change_names <- function(data, var = "Outcome", new = c("OS probability", "DFS probability")){

    levels(data[[var]]) <- new
 return(data)
  }

#


KM_GCjohn <- subset(kmdata, Estimate == "Orig. IPD" | Estimate == "GC all mom.")

mp <- aes( x= time, y= km)
cis <- geom_ribbon(aes( ymin=km_l, ymax=km_u, fill = Treatment ), alpha=0.5)
panel <- facet_grid( Outcome ~ Estimate)
textsize <- theme(axis.text = element_text(size=18), axis.title = element_text(size=20), strip.text.x = element_text(size = 19), strip.text.y = element_text(size = 19), legend.title=element_blank(), legend.text=element_text(size=18), legend.position="bottom" )

#
 png(file.path(outfolder, "km_plot_meta_CIs.png"), width = 1000, height = 800)

 ggplot( change_names(KM_MA), mp ) + cis + panel + xl + yl + xlims + theme_bw() + textsize 

 dev.off()

#
 png(file.path(outfolder, "km_plot_GCgamma_CIs.png"), width = 1000, height = 800)

 ggplot( change_names(KM_GC_gam), mp ) + cis + panel + xl + yl + xlims + theme_bw() + textsize 

dev.off()

                                        #
 png(file.path(outfolder, "km_plot_GCjohn_CIs.png"), width = 1000, height = 800)

 ggplot( change_names(KM_GCjohn), mp ) + cis + panel + xl + yl + xlims + theme_bw() + textsize 

 dev.off()

                                        #
 png(file.path(outfolder, "km_plot_synpd_CIs.png"), width = 1000, height = 800)

 ggplot( change_names(KM_CART), mp ) + cis + panel + xl + yl + xlims + theme_bw() + textsize 

 dev.off()

                                        #
 png(file.path(outfolder, "km_plot_vae_CIs.png"), width = 1000, height = 800)

 ggplot( change_names(KM_VAE), mp ) + cis + panel + xl + yl + xlims + theme_bw() + textsize 

 dev.off()


            #
 png(file.path(outfolder, "km_plot_bootDist_CIs.png"), width = 1000, height = 800)

 ggplot( change_names(KM_DISTBOOT), mp ) + cis + panel + xl + yl + xlims + theme_bw() + textsize 

 dev.off()
