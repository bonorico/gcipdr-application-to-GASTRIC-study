## other_pseudodata_methods.R contains R commands to execute 'gcipdr application to GASTRIC group meta-analysis' (from homonymous repository). 
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



############## SYNTHPOP COMPARISON (SUPPLEMENTARY MATERIAL)

# each center sends synth data: default method (CART gen, simple sampling)

system.time(
    regio.synpopdat <- lapply(trialnames,
                              function(tn)
                                  syn(gastadj[gastadj$trialID == tn, -1],
                                      seed =  0125 + as.integer(tn),
                                      m = H
                                      )
                              )
)

names(regio.synpopdat) <- trialnames
      
  
                                        # pool centrally

pool.synpdat <- lapply(1:H,
                       function(h)
                           do.call(
                               "rbind",
                               lapply(trialnames,
                                      function(tn)
                                          data.frame(regio.synpopdat[[tn]]$syn[[h]],
                                                     trialID = tn)

                                      )
                           )
                       )



############## VAE COMPARISON ###############


system.time(
    regio.vae <- lapply(trialnames,
                        function(tn)
                        {
                            set.seed(472 + as.integer(tn) )
                            print(
                                system.time(
                                    vaeout <- VAE(
                                        gastadj[gastadj$trialID == tn, -1],
                                        H,
                                        epochs = 50,
                                        repeat_data = 100
                                    )

                                )
                            )
                            return(vaeout)

                        }
                        )

)

names(regio.vae) <- trialnames


# vaeres <- VAE(gastadj[gastadj$trialID == "1", -1], 100, repeat_data = 100, epochs = 50 )
## TODO: check vae again ... line above yields different output than below (where ?) .... why ? BECAUSE set.seed does not work within lapply

pool.vaedat <- lapply(1:H,
                      function(h)
                          do.call(
                              "rbind",
                              lapply(trialnames,
                                     function(tn)
                                         data.frame(regio.vae[[tn]][[h]],
                                                    trialID = tn)
                                     )
                          )
                      )


############## IPD BOOTSTRAP COMPARISON ###############

set.seed(589)
bootres <- IPDbootDat(gastadj, H, distributed = FALSE, site = "trialID")


############## DISTRIBUTED IPD BOOTSTRAP COMPARISON ###############

set.seed(588)
bootresDist <- IPDbootDat(gastadj, H, distributed = TRUE, site = "trialID")



print("New objects: pool.synpdat, pool.vaedat, bootres, bootresDist")
