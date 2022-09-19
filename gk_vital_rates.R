#=====================
# Killifish-Guppy vital rates 
#    MAIN SCRIPT
#=====================
# Uploads data, fits Bayesian vital rate models simultaneously using NIMBLE
# model diagnostics, model selection, etc...
#=====================
#Load necessary libraries and functions
library(Rcpp)
library(devtools)
require(nimble)
library(plyr)
library(MCMCvis)

source('./upload_data.R')

#===============================================================================================
#Model with growth = log(z'/z), phis, rhos, interaction with body size, single-pop channels only 
#===============================================================================================

# MCMC parameters 
ni <- 1200000 # total iterations per chain
nb <-  200000 # burn-in period
nt <- 1000 # thinning interval (to save memory)
nc <- 4 # number of chains to run

params <- c('gKG.beta','gKO.beta','gLP.beta','gHP.beta',
            'fKG.beta','fKO.beta','fLP.beta','fHP.beta',
            'philength.KG', 'philength.KO','philength.LP','philength.HP',
            'rholength.KG', 'rholength.KO','rholength.LP','rholength.HP',
            'g.sigma', 'g.sigma.int.c')


gdata <- upload_data() # upload data
source('./gk_model3.R') # upload model code

# build the model
Rmodel <- nimbleModel(code1, constants, data, inits)
conf <- configureMCMC(Rmodel, monitors = params, enableWAIC = TRUE)
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel,showCompilerOutput = TRUE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel,showCompilerOutput = TRUE)

# fit the model
set.seed(0)
samples_model3 <- runMCMC(Cmcmc, niter = ni, nburnin = nb, thin = nt,
                                nchains = nc, WAIC = TRUE)


save(samples_model3, file="samples_model3.Rda")


# 
# # plots of posteriors and trace to check convergence etc
# MCMCtrace(samples_model3$samples, 
#           params = params, pdf=F, ind = T, Rhat = T, n.eff = T)
# 
# # summary of parameter estimates
# samplesSummary(samples_model2$samples[[1]])
# 
# # get WAIC
# samples_model4$WAIC
# 
# save(samples_model3, file="samples_model3.Rda")
# 
# 
# 
# 
# ni <- 120000 # total iterations per chain
# nb <-  20000 # burn-in period
# nt <- 100 # thinning interval (to save memory)
# nc <- 2 # number of chains to run
# 
# params <- c('gKG.beta','gKO.beta','gLP.beta','gHP.beta',
#             'fKG.beta','fKO.beta','fLP.beta','fHP.beta',
#             'philength.KG', 'philength.KO','philength.LP','philength.HP',
#             'rholength.KG', 'rholength.KO','rholength.LP','rholength.HP',
#             'g.sigma', 'g.sigma.int.c')
# 
# 
# gdata <- upload_data() # upload data
# source('./gk_model2.R') # upload model code
# 
# # build the model
# Rmodel <- nimbleModel(code1, constants, data, inits)
# conf <- configureMCMC(Rmodel, monitors = params, enableWAIC = TRUE)
# Rmcmc <- buildMCMC(conf)
# Cmodel <- compileNimble(Rmodel,showCompilerOutput = TRUE)
# Cmcmc <- compileNimble(Rmcmc, project = Rmodel,showCompilerOutput = TRUE)
# 
# # fit the model
# set.seed(0)
# samples_model1 <- runMCMC(Cmcmc, niter = ni, nburnin = nb, thin = nt,
#                           nchains = nc, WAIC = TRUE)
# 
# 
# # plots of posteriors and trace to check convergence etc
# MCMCtrace(samples_model1$samples, 
#           params = params, pdf=F, ind = T, Rhat = T, n.eff = T)
# 
# # summary of parameter estimates
# samplesSummary(samples_model1$samples[[1]])
# 
# # get WAIC
# samples_model4$WAIC
# 
# save(samples_model1, file="samples_model1.Rda")
# 
# 
# library(glmmTMB)
# 
# test <- gdata
# test$ln_init_sl <- log(test$init.sl/6.5)
# test$dens <- test$total.density/1.4
# 
# mod_gKG <- glmmTMB(growth ~ 1 + ln_init_sl + dens + dens:ln_init_sl +
#                      (1|channel.num),
#                    data=test[KG==1,])
# 
# mod_gKO <- glmmTMB(growth ~ 1 + ln_init_sl + dens + dens:ln_init_sl +
#                      (1|channel.num),
#                    data=test[KO==1,])
# 
# mod_gLP <- glmmTMB(growth ~ 1 + ln_init_sl + dens + dens:ln_init_sl +
#                      (1|channel.num),
#                    data=test[LP==1,])
# 
# mod_gHP <- glmmTMB(growth ~ 1 + ln_init_sl + dens + dens:ln_init_sl +
#                      (1|channel.num),
#                    data=test[HP==1,])
# 
# #fecundity
# mod_fKG <- glmmTMB(num.embryos ~ 1 + ln_init_sl + dens + dens:ln_init_sl, 
#                    family="poisson",
#                    data=test[KG==1,])
# 
# mod_fKO <- glmmTMB(num.embryos ~ 1 + ln_init_sl + dens + dens:ln_init_sl,
#                    family="poisson",
#                    data=test[KO==1,])
# 
# mod_fLP <- glmmTMB(num.embryos ~ 1 + ln_init_sl + dens + dens:ln_init_sl,
#                    family="poisson",
#                    data=test[LP==1,])
# 
# mod_fHP <- glmmTMB(num.embryos ~ 1 + ln_init_sl + dens + dens:ln_init_sl,
#                    family="poisson",
#                    data=test[HP==1,])
