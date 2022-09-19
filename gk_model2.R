
#These functions are required in the likelihood function. They create the interaction surface.

length.cen <- 6.5

# # KG kilifish focal with KG killifish competitors
# myFunctionKGKG <- nimbleFunction(
#   run = function(H = double(0), KG = double(0), philength.KG = double(0), rholength.KG = double(0), density.KG = double(0), le.comp.KG = double(1), le.focal.KGKG = double(1)) {
#     numNumber <- KG * (sum( ((exp(philength.KG * le.comp.KG - philength.KG * le.focal.KGKG) ) *
#                                exp((-(rholength.KG * le.comp.KG - rholength.KG * le.focal.KGKG)^2) / (2*H^2) ))[1:density.KG]
#     ) ) / 1.4
#     
#     returnType(double(0))
#     return(numNumber)
#   }
# )
# 
# # KO kilifish with KO killifish
# myFunctionKOKO <- nimbleFunction(
#   run = function(H = double(0), KO = double(0), philength.KO = double(0), rholength.KO = double(0), density.KO = double(0), le.comp.KO = double(1), le.focal.KOKO = double(1)) {
#     numNumber <- KO * (sum( ((exp(philength.KO * le.comp.KO - philength.KO * le.focal.KOKO) ) *
#                                exp((-(rholength.KO * le.comp.KO - rholength.KO * le.focal.KOKO)^2) / (2*H^2) ))[1:density.KO]
#     ) ) / 1.4
#     
#     returnType(double(0))
#     return(numNumber)
#   }
# )
# 
# 
# # LP guppies focal with LP guppy competitors
# myFunctionLPLP <- nimbleFunction(
#   run = function(H = double(0), LP = double(0), philength.LP = double(0), rholength.LP = double(0), density.LP = double(0), le.comp.LP = double(1), le.focal.LPLP = double(1)) {
#     numNumber <- LP * (sum( ((exp(philength.LP * le.comp.LP - philength.LP * le.focal.LPLP) ) *
#                                exp((-(rholength.LP * le.comp.LP - rholength.LP * le.focal.LPLP)^2) / (2*H^2) ))[1:density.LP]
#     ) ) / 1.4
#     
#     returnType(double(0))
#     return(numNumber)
#   }
# )
# 
# # HP guppies with HP guppies
# myFunctionHPHP <- nimbleFunction(
#   run = function(H = double(0), HP = double(0), philength.HP = double(0), rholength.HP = double(0), density.HP = double(0), le.comp.HP = double(1), le.focal.HPHP = double(1)) {
#     numNumber <- HP * (sum( ((exp(philength.HP * le.comp.HP - philength.HP * le.focal.HPHP) ) *
#                                exp((-(rholength.HP * le.comp.HP - rholength.HP * le.focal.HPHP)^2) / (2*H^2) ))[1:density.HP]
#     ) ) / 1.4
#     
#     returnType(double(0))
#     return(numNumber)
#   }
# )
# 
# 
# # KG killifish focal with LP guppy competitor
# myFunctionKGLP <- nimbleFunction(
#   run = function(H = double(0), KG = double(0), philength.KG = double(0), rholength.KG = double(0),
#                  philength.LP = double(0), rholength.LP = double(0),
#                  density.KG = double(0), le.comp.KG = double(1), le.comp.LP = double(1),le.focal.KGKG = double(1),
#                  iota.KiGu = double(0), eta.KiGu = double(0))
#     {
#     numNumber <- KG * (sum( ((exp(iota.KiGu + philength.LP * le.comp.LP - philength.KG * le.focal.KGKG) ) *
#                                exp((-(eta.KiGu + rholength.LP * le.comp.LP - rholength.KG * le.focal.KGKG)^2) / (2*H^2) ))[1:density.KG]
#     ) ) / 1.4
#     
#     returnType(double(0))
#     return(numNumber)
#   }
# )
# # KG killifish focal with HP guppy competitor 
# myFunctionKGHP <- nimbleFunction(
#   run = function(H = double(0), KG = double(0), philength.KG = double(0), rholength.KG = double(0),
#                  philength.HP = double(0), rholength.HP = double(0),
#                  density.KG = double(0), le.comp.KG = double(1), le.comp.HP = double(1),le.focal.KGKG = double(1),
#                  iota.KiGu = double(0), eta.KiGu = double(0))
#   {
#     numNumber <- KG * (sum( ((exp(iota.KiGu + philength.HP * le.comp.HP - philength.KG * le.focal.KGKG) ) *
#                                exp((-(eta.KiGu + rholength.HP * le.comp.HP - rholength.KG * le.focal.KGKG)^2) / (2*H^2) ))[1:density.KG]
#     ) ) / 1.4
#     
#     returnType(double(0))
#     return(numNumber)
#   }
# )
# 
# # KO killifish focal with LP guppy competitor
# myFunctionKOLP <- nimbleFunction(
#   run = function(H = double(0), KO = double(0), philength.KO = double(0), rholength.KO = double(0),
#                  philength.LP = double(0), rholength.LP = double(0),
#                  density.KO = double(0), le.comp.KO = double(1), le.comp.LP = double(1),le.focal.KOKO = double(1),
#                  iota.KiGu = double(0), eta.KiGu = double(0))
#   {
#     numNumber <- KO * (sum( ((exp(iota.KiGu + philength.LP * le.comp.LP - philength.KO * le.focal.KOKO) ) *
#                                exp((-(eta.KiGu + rholength.LP * le.comp.LP - rholength.KO * le.focal.KOKO)^2) / (2*H^2) ))[1:density.KO]
#     ) ) / 1.4
#     
#     returnType(double(0))
#     return(numNumber)
#   }
# )
# # KO killifish focal with HP guppy competitor 
# myFunctionKOHP <- nimbleFunction(
#   run = function(H = double(0), KO = double(0), philength.KO = double(0), rholength.KO = double(0),
#                  philength.HP = double(0), rholength.HP = double(0),
#                  density.KO = double(0), le.comp.KO = double(1), le.comp.HP = double(1),le.focal.KOKO = double(1),
#                  iota.KiGu = double(0), eta.KiGu = double(0))
#   {
#     numNumber <- KO * (sum( ((exp(iota.KiGu + philength.HP * le.comp.HP - philength.KO * le.focal.KOKO) ) *
#                                exp((-(eta.KiGu + rholength.HP * le.comp.HP - rholength.KO * le.focal.KOKO)^2) / (2*H^2) ))[1:density.KO]
#     ) ) / 1.4
#     
#     returnType(double(0))
#     return(numNumber)
#   }
# )
# 
# # LP guppy focal with KG killifish competitor
# myFunctionLPKG <- nimbleFunction(
#   run = function(H = double(0), LP = double(0), philength.LP = double(0), rholength.LP = double(0),
#                  philength.KG = double(0), rholength.KG = double(0),
#                  density.LP = double(0), le.comp.LP = double(1), le.comp.KG = double(1),le.focal.LPLP = double(1),
#                  iota.KiGu = double(0), eta.KiGu = double(0))
#   {
#     numNumber <- LP * (sum( ((exp(-iota.KiGu + philength.KG * le.comp.KG - philength.LP * le.focal.LPLP) ) *
#                                exp((-(-eta.KiGu + rholength.KG * le.comp.KG - rholength.LP * le.focal.LPLP)^2) / (2*H^2) ))[1:density.LP]
#     ) ) / 1.4
#     
#     returnType(double(0))
#     return(numNumber)
#   }
# )
# 
# # LP guppy focal with KO killifish competitor
# myFunctionLPKO <- nimbleFunction(
#   run = function(H = double(0), LP = double(0), philength.LP = double(0), rholength.LP = double(0),
#                  philength.KO = double(0), rholength.KO = double(0),
#                  density.LP = double(0), le.comp.LP = double(1), le.comp.KO = double(1),le.focal.LPLP = double(1),
#                  iota.KiGu = double(0), eta.KiGu = double(0))
#   {
#     numNumber <- LP * (sum( ((exp(-iota.KiGu + philength.KO * le.comp.KO - philength.LP * le.focal.LPLP) ) *
#                                exp((-(-eta.KiGu + rholength.KO * le.comp.KO - rholength.LP * le.focal.LPLP)^2) / (2*H^2) ))[1:density.LP]
#     ) ) / 1.4
#     
#     returnType(double(0))
#     return(numNumber)
#   }
# )
# 
# # HP guppy focal with KG killifish competitor
# myFunctionHPKG <- nimbleFunction(
#   run = function(H = double(0), HP = double(0), philength.HP = double(0), rholength.HP = double(0),
#                  philength.KG = double(0), rholength.KG = double(0),
#                  density.HP = double(0), le.comp.HP = double(1), le.comp.KG = double(1),le.focal.HPHP = double(1),
#                  iota.KiGu = double(0), eta.KiGu = double(0))
#   {
#     numNumber <- HP * (sum( ((exp(-iota.KiGu + philength.KG * le.comp.KG - philength.HP * le.focal.HPHP) ) *
#                                exp((-(-eta.KiGu + rholength.KG * le.comp.KG - rholength.HP * le.focal.HPHP)^2) / (2*H^2) ))[1:density.HP]
#     ) ) / 1.4
#     
#     returnType(double(0))
#     return(numNumber)
#   }
# )
# 
# # HP guppy focal with KO killifish competitor
# myFunctionHPKO <- nimbleFunction(
#   run = function(H = double(0), HP = double(0), philength.HP = double(0), rholength.HP = double(0),
#                  philength.KO = double(0), rholength.KO = double(0),
#                  density.HP = double(0), le.comp.HP = double(1), le.comp.KO = double(1),le.focal.HPHP = double(1),
#                  iota.KiGu = double(0), eta.KiGu = double(0))
#   {
#     numNumber <- HP * (sum( ((exp(-iota.KiGu + philength.KO * le.comp.KO - philength.HP * le.focal.HPHP) ) *
#                                exp((-(-eta.KiGu + rholength.KO * le.comp.KO - rholength.HP * le.focal.HPHP)^2) / (2*H^2) ))[1:density.HP]
#     ) ) / 1.4
#     
#     returnType(double(0))
#     return(numNumber)
#   }
# )

# GENERAL INTERACTION SURFACCE FUNCTION

competitionSurface <- nimbleFunction(
  run = function(H = double(0), population = double(0), philength.focal = double(0), rholength.focal = double(0),
                 philength.comp = double(0), rholength.comp = double(0),
                 density.focal = double(0), le.comp = double(1), le.focal = double(1),
                 iota = double(0), eta = double(0))
  {
    # numNumber <- pop * (sum( ((exp(iota + philength.comp * le.comp - philength.focal * le.focal) ) *
    #                            exp((-(eta + rholength.comp * le.comp - rholength.focal * le.focal)^2) / (2*H^2) ))[1:density.focal]
    # ) ) / 1.4
    
    numNumber <- population  * (sum( ((exp(0 + 0 * le.comp - 0 * le.focal) ) *
                                exp((-(0 + 0 * le.comp - 0 * le.focal)^2) / (2*1^2) ))[1:density.focal]
    ) ) / 1.4
    
    returnType(double(0))
    return(numNumber)
  }
)

#This is the Nimble model
code1 <- nimbleCode({
  #=====================================================================================
  #Priors for betas
  #=====================================================================================
  for(j in 1:4) {
    gKG.beta[j] ~ dnorm(0, sd = 100)
    gKO.beta[j] ~ dnorm(0, sd = 100)
    gLP.beta[j] ~ dnorm(0, sd = 100)
    gHP.beta[j] ~ dnorm(0, sd = 100)
  }
  for(j in 1:4) {
    fKG.beta[j] ~ dnorm(0, sd = 100)
    fKO.beta[j] ~ dnorm(0, sd = 100)
    fLP.beta[j] ~ dnorm(0, sd = 100)
    fHP.beta[j] ~ dnorm(0, sd = 100)
    
  }
  
  
  
  #=====================================================================================
  #Priors for interaction surface
  #=====================================================================================
  
  philength.KG ~ dnorm(0.067, 0.010) #From Analysis with Interaction
  philength.KO ~ dnorm(0.020, 0.017) #From Analysis with Interaction
  philength.LP ~ dnorm(0.05, 0.009) #From Analysis with Interaction
  philength.HP ~ dnorm(0.02, 0.010) #From Analysis with Interaction
  
  rholength.KG ~ dnorm(0.03, 0.005) #From Analysis with Interaction
  rholength.KO ~ dnorm(0.03, 0.004) #From Analysis with Interaction
  rholength.LP ~ dnorm(-0.01, 0.017) #From Analysis with Interaction
  rholength.HP ~ dnorm(0.10, 0.021) #From Analysis with Interaction
  
  
  g.sigma ~ dunif(0, 100)
  
  #=====================================================================================
  #Priors for random effects
  #=====================================================================================
  
  ## random intercepts
  # random intercept of channel on growth
  g.sigma.int.c ~ dunif(0, 100)
  g.tau.int.c <- 1/(g.sigma.int.c*g.sigma.int.c)
  for(ii in 1:n_chan) {
    g.alpha.channel[ii] ~ dnorm(0, g.tau.int.c)
  }
  
  #=====================================================================================
  #Likelihood
  #=====================================================================================
  
  for(i in 1:n_obs) {
    
    # define the interaction surface for interactions between focal and competitors:
    
    intsurf.KGKG[i] <-  competitionSurface(H=H, population=KG[i], philength.focal=philength.KG, rholength.focal=rholength.KG, 
                                           philength.comp=philength.KG, rholength.comp=rholength.KG,
                                           density.focal=density.KG[i], le.comp=le.comp.KG[, chan.num[i]], le.focal=le.focal.KG[i,],
                                           iota=0, eta=0)
    
    intsurf.KOKO[i] <-  competitionSurface(H=H, population=KO[i], philength.focal=philength.KO, rholength.focal=rholength.KO, 
                                           philength.comp=philength.KO, rholength.comp=rholength.KO,
                                           density.focal=density.KO[i], le.comp=le.comp.KO[, chan.num[i]], le.focal=le.focal.KO[i,],
                                           iota=0, eta=0)
    
    intsurf.LPLP[i] <-  competitionSurface(H=H, population=LP[i], philength.focal=philength.LP, rholength.focal=rholength.LP,
                                           philength.comp=philength.LP, rholength.comp=rholength.LP, 
                                           density.focal=density.LP[i], le.comp=le.comp.LP[, chan.num[i]], le.focal=le.focal.LP[i,],
                                           iota=0, eta=0)

    intsurf.HPHP[i] <-  competitionSurface(H=H, population=HP[i], philength.focal=philength.HP, rholength.focal=rholength.HP, 
                                           philength.comp=philength.HP, rholength.comp=rholength.HP, 
                                           density.focal=density.HP[i], le.comp=le.comp.HP[, chan.num[i]], le.focal=le.focal.HP[i,],
                                           iota=0, eta=0)
    # 
    # intsurf.KGLP[i] <-  competitionSurface(H=H, pop=KG[i], philength.focal=philength.KG, rholength.focal=rholength.KG, 
    #                                        philength.comp=philength.LP, rholength.comp=rholength.LP,
    #                                        density.focal=density.KG[i], le.comp=le.comp.LP[, chan.num[i]], le.focal=le.focal.LP[i,],
    #                                        iota=iota.KiGu, eta=eta.KiGu)
    # 
    # intsurf.KGHP[i] <-  competitionSurface(H=H, pop=KG[i], philength.focal=philength.KG, rholength.focal=rholength.KG, 
    #                                        philength.comp=philength.HP, rholength.comp=rholength.HP,
    #                                        density.focal=density.KG[i], le.comp=le.comp.HP[, chan.num[i]], le.focal=le.focal.HP[i,],
    #                                        iota=iota.KiGu, eta=eta.KiGu)
    # 
    # intsurf.KOLP[i] <-  competitionSurface(H=H, pop=KO[i], philength.focal=philength.KO, rholength.focal=rholength.KO, 
    #                                        philength.comp=philength.LP, rholength.comp=rholength.LP,
    #                                        density.focal=density.KO[i], le.comp=le.comp.LP[, chan.num[i]], le.focal=le.focal.LP[i,],
    #                                        iota=iota.KiGu, eta=eta.KiGu)
    # 
    # intsurf.KOHP[i] <-  competitionSurface(H=H, pop=KO[i], philength.focal=philength.KO, rholength.focal=rholength.KO, 
    #                                        philength.comp=philength.HP, rholength.comp=rholength.HP,
    #                                        density.focal=density.KO[i], le.comp=le.comp.HP[, chan.num[i]], le.focal=le.focal.HP[i,],
    #                                        iota=iota.KiGu, eta=eta.KiGu)
    # 
    # intsurf.LPKG[i] <-  competitionSurface(H=H, pop=LP[i], philength.focal=philength.LP, rholength.focal=rholength.LP, 
    #                                        philength.comp=philength.KG, rholength.comp=rholength.KG,
    #                                        density.focal=density.LP[i], le.comp=le.comp.KG[, chan.num[i]], le.focal=le.focal.KG[i,],
    #                                        iota=-iota.KiGu, eta=-eta.KiGu)
    # 
    # intsurf.LPKO[i] <-  competitionSurface(H=H, pop=LP[i], philength.focal=philength.LP, rholength.focal=rholength.LP, 
    #                                        philength.comp=philength.KO, rholength.comp=rholength.KO,
    #                                        density.focal=density.LP[i], le.comp=le.comp.KO[, chan.num[i]], le.focal=le.focal.KO[i,],
    #                                        iota=-iota.KiGu, eta=-eta.KiGu)
    # 
    # intsurf.HPKG[i] <-  competitionSurface(H=H, pop=HP[i], philength.focal=philength.HP, rholength.focal=rholength.HP, 
    #                                        philength.comp=philength.KG, rholength.comp=rholength.KG,
    #                                        density.focal=density.HP[i], le.comp=le.comp.KG[, chan.num[i]], le.focal=le.focal.KG[i,],
    #                                        iota=-iota.KiGu, eta=-eta.KiGu)
    # 
    # intsurf.HPKO[i] <-  competitionSurface(H=H, pop=HP[i], philength.focal=philength.HP, rholength.focal=rholength.HP, 
    #                                        philength.comp=philength.KO, rholength.comp=rholength.KO,
    #                                        density.focal=density.HP[i], le.comp=le.comp.KO[, chan.num[i]], le.focal=le.focal.KO[i,],
    #                                        iota=-iota.KiGu, eta=-eta.KiGu)
    
    
    
    g.mu[i] <- KO[i] * (gKO.beta[1] + gKO.beta[2]*initial.length[i] + (gKO.beta[3] + gKO.beta[4]*initial.length[i])*intsurf.KOKO[i]) +
                                                          # KO.LP[i] * (gKO.beta[3] + gKO.beta[4]*initial.length[i])*intsurf.KOLP[i] +
                                                          # KO.HP[i] * (gKO.beta[3] + gKO.beta[4]*initial.length[i])*intsurf.KOHP[i]) +
                          
               KG[i] * (gKG.beta[1] + gKG.beta[2]*initial.length[i] + (gKG.beta[3] + gKG.beta[4]*initial.length[i])*intsurf.KGKG[i]) +
                                                          # KG.LP[i] * (gKG.beta[3] + gKG.beta[4]*initial.length[i])*intsurf.KGLP[i] +
                                                          # KG.HP[i] * (gKG.beta[3] + gKG.beta[4]*initial.length[i])*intsurf.KGHP[i]) +
      
               LP[i] * (gLP.beta[1] + gLP.beta[2]*initial.length[i] + (gLP.beta[3] + gLP.beta[4]*initial.length[i])*intsurf.LPLP[i]) +
                                                          # KG.LP[i] * (gLP.beta[3] + gLP.beta[4]*initial.length[i])*intsurf.LPKG[i] +
                                                          # KO.LP[i] * (gLP.beta[3] + gLP.beta[4]*initial.length[i])*intsurf.LPKO[i]) +
      
               HP[i] * (gHP.beta[1] + gHP.beta[2]*initial.length[i] + (gHP.beta[3] + gHP.beta[4]*initial.length[i])*intsurf.HPHP[i]) +
                                                          # KG.HP[i] * (gHP.beta[3] + gHP.beta[4]*initial.length[i])*intsurf.HPKG[i] +
                                                           #KO.HP[i] * (gHP.beta[3] + gHP.beta[4]*initial.length[i])*intsurf.HPKO[i]) +
      
               g.alpha.channel[chan.num[i]] 
    
    growth[i] ~ dnorm(g.mu[i],sd=g.sigma)
    
    
    log(f.mu[i]) <-  KO[i] * (fKO.beta[1] + fKO.beta[2]*initial.length[i] + (fKO.beta[3] + fKO.beta[4]*initial.length[i])*intsurf.KOKO[i]) +
                                                               #  KO.LP[i] * (fKO.beta[3] + fKO.beta[4]*initial.length[i])*intsurf.KOLP[i] +
                                                              #   KO.HP[i] * (fKO.beta[3] + fKO.beta[4]*initial.length[i])*intsurf.KOHP[i]) +
      
                     KG[i] * (fKG.beta[1] + fKG.beta[2]*initial.length[i] + (fKG.beta[3] + fKG.beta[4]*initial.length[i])*intsurf.KGKG[i]) +
                                                               #  KG.LP[i] * (fKG.beta[3] + fKG.beta[4]*initial.length[i])*intsurf.KGLP[i] +
                                                                # KG.HP[i] * (fKG.beta[3] + fKG.beta[4]*initial.length[i])*intsurf.KGHP[i]) +
                     
                     LP[i] * (fLP.beta[1] + fLP.beta[2]*initial.length[i] + (fLP.beta[3] + fLP.beta[4]*initial.length[i])*intsurf.LPLP[i]) +
                                                                # KG.LP[i] * (fLP.beta[3] + fLP.beta[4]*initial.length[i])*intsurf.LPKG[i] +
                                                                 #KO.LP[i] * (fLP.beta[3] + fLP.beta[4]*initial.length[i])*intsurf.LPKO[i]) +
                     
                     HP[i] * (fHP.beta[1] + fHP.beta[2]*initial.length[i] + (fHP.beta[3] + fHP.beta[4]*initial.length[i])*intsurf.HPHP[i])# +
                                                                # KG.HP[i] * (fHP.beta[3] + fHP.beta[4]*initial.length[i])*intsurf.HPKG[i] +
                                                                # KO.HP[i] * (fHP.beta[3] + fHP.beta[4]*initial.length[i])*intsurf.HPKO[i])
    
    
    fec[i] ~ dpois(f.mu[i])
    
    
    # d.mu[i] <-  KO[i] * (d.beta[1] ) +
    #    KG[i] * (d.beta[4] ) 
    
    
    # off.length[i] ~ dnorm(d.mu[i],sd=d.sigma)
    
    
  }
})


##########################################################
#These are required to make the interaction surface
##########################################################
# pre-allocate a list and fill it with a loop

# specify focal traits

le.focal.KG <- array(0,c(dim(gdata)[1],24))
for (k in 1:dim(gdata)[1]){
  i <- gdata[,channel.num[k]]
  if (gdata[, KG[k]==1]){
    le.focal.KG[k,c(1:gdata[channel.num==i & KG==1, length(init.sl)])] <- gdata[,init.sl[k]-length.cen]
  }
}

le.focal.KO <- array(0,c(dim(gdata)[1],24))
for (k in 1:dim(gdata)[1]){
  i <- gdata[,channel.num[k]]
  if (gdata[, KO[k]==1]){
    le.focal.KO[k,c(1:gdata[channel.num==i & KO==1, length(init.sl)])] <- gdata[,init.sl[k]-length.cen]
  }
}

le.focal.LP <- array(0,c(dim(gdata)[1],24))
for (k in 1:dim(gdata)[1]){
  i <- gdata[,channel.num[k]]
  if (gdata[, LP[k]==1]){
    le.focal.LP[k,c(1:gdata[channel.num==i & LP==1, length(init.sl)])] <- gdata[,init.sl[k]-length.cen]
  }
}

le.focal.HP <- array(0,c(dim(gdata)[1],24))
for (k in 1:dim(gdata)[1]){
  i <- gdata[,channel.num[k]]
  if (gdata[, HP[k]==1]){
    le.focal.HP[k,c(1:gdata[channel.num==i & HP==1, length(init.sl)])] <- gdata[,init.sl[k]-length.cen]
  }
}

# specify competitor traits

le.comp.KG <- array(0,c(24,length(unique(gdata$channel.num))))
for (i in 1:max(gdata$channel.num)){
  for (j in 1:gdata[channel.num==i,length(init.sl)]){
    le.comp.KG[j,i] <- gdata[channel.num==i & KG==1, init.sl[j]-length.cen] 
  }
}

le.comp.KO <- array(0,c(24,length(unique(gdata$channel.num))))
for (i in 1:max(gdata$channel.num)){
  for (j in 1:gdata[channel.num==i,length(init.sl)]){
    le.comp.KO[j,i] <- gdata[channel.num==i & KO==1, init.sl[j]-length.cen] 
  }
}

le.comp.LP <- array(0,c(24,length(unique(gdata$channel.num))))
for (i in 1:max(gdata$channel.num)){
  for (j in 1:gdata[channel.num==i,length(init.sl)]){
    le.comp.LP[j,i] <- gdata[channel.num==i & LP==1, init.sl[j]-length.cen] 
  }
}

le.comp.HP <- array(0,c(24,length(unique(gdata$channel.num))))
for (i in 1:max(gdata$channel.num)){
  for (j in 1:gdata[channel.num==i,length(init.sl)]){
    le.comp.HP[j,i] <- gdata[channel.num==i & HP==1, init.sl[j]-length.cen] 
  }
}

le.comp.KG[is.na(le.comp.KG)==T] <- 0
le.comp.KO[is.na(le.comp.KO)==T] <- 0
le.comp.LP[is.na(le.comp.LP)==T] <- 0
le.comp.HP[is.na(le.comp.HP)==T] <- 0







#growth = gdata[,growth] 
#fec = gdata[,num.embryos]

#n_block = gdata[,max(as.numeric(block))]
n_chan = gdata[,max(channel.num)]
# n_obs = dim(gdata)[1]
# chan.num = gdata[,channel.num] 
# block.num = gdata[,block] 
# #exper.num=as.numeric(gdata[,'exper'])
# initial.length = gdata[, init.sl - length.cen] 
# #initial.length.squ = (as.numeric(gdata[,'init.sl'] - length.cen))^2
# #stage = as.numeric(gdata[,'stage'])
# le.comp.KG = le.comp.KG
# le.comp.KO = le.comp.KO
# le.comp.LP = le.comp.LP
# le.comp.HP = le.comp.HP
# le.focal.KGKG = le.focal.KGKG
# le.focal.KOKO = le.focal.KOKO
# le.focal.LPLP = le.focal.LPLP
# le.focal.HPHP = le.focal.HPHP
# # le.focal.KOKG = le.focal.KOKG
# # le.focal.KGKO = le.focal.KGKO
# KG = gdata[,KG] 
# KO = gdata[,KO]
# LP = gdata[,LP] 
# HP = gdata[,HP]
# # density=as.numeric(gdata[,'density'])
# density.KG = gdata[,density.KG]
# #density.KG[which(density.KG==0)] <- 1
# density.KO = gdata[,density.KO]
# #density.KO[which(density.KO==0)] <- 1
# density.LP = gdata[,density.LP]
# #density.LP[which(density.LP==0)] <- 1
# density.HP = gdata[,density.HP]
# #density.HP[which(density.HP==0)] <- 1

# density = density.KG + density.KO
#both.pheno = gdata[,both.pheno]


# specify the predictors
constants <- list(
  chan.num = gdata[,channel.num],
  block.num = gdata[,block],
  n_obs = dim(gdata)[1],
  initial.length = gdata[, log(init.sl / length.cen)],
  # initial.length.fec = initial.length.fec,
  #  initial.length.squ = initial.length.squ,
  # stage = stage,
  KO = gdata[,KO],
  KG = gdata[,KG],
  LP = gdata[,LP],
  HP = gdata[,HP],
  # both.pheno = both.pheno,
  density.KG = gdata[,density.KG],
  density.KO = gdata[,density.KO],
  density.LP = gdata[,density.LP],
  density.HP = gdata[,density.HP],
  le.comp.KG = le.comp.KG,
  le.comp.KO = le.comp.KO,
  le.comp.LP = le.comp.LP,
  le.comp.HP = le.comp.HP,
  le.focal.KG = le.focal.KG,
  le.focal.KO = le.focal.KO,
  le.focal.LP = le.focal.LP,
  le.focal.HP = le.focal.HP
  #le.focal.KOKG = le.focal.KOKG,
  #le.focal.KGKO = le.focal.KGKO,
)

# specify the observations
data <- list(
  growth = gdata[,growth],
  fec = gdata[,num.embryos]
)

# specify initial values for all model parameters
inits <- list(
  gKO.beta = rep(0, 4),
  gKG.beta = rep(0, 4),
  gLP.beta = rep(0, 4),
  gHP.beta = rep(0, 4),
  fKO.beta = rep(0, 4),
  fKG.beta = rep(0, 4),
  fLP.beta = rep(0, 4),
  fHP.beta = rep(0, 4),
  
  philength.KG = 0,
  philength.KO = 0,
  philength.LP = 0,
  philength.HP = 0,
  rholength.KG = 0,
  rholength.KO = 0,
  rholength.LP = 0,
  rholength.HP = 0,
  H=1,
  g.sigma = 1,
  g.sigma.int.c = 1,
  g.alpha.channel = rep(0, n_chan)
  
)

