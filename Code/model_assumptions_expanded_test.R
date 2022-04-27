# model assumptions test
setwd("~/+My_Documents/MRSA risk group model/MRSA Risk Group Model Code")
.libPaths("C:/Users/nrt8/Documents/R/win-library/4.0")
source("model_functions.R")
require(deSolve)
require(rootSolve)
require(tidyverse)
require(cowplot)

# set plot font size for everything
theme_set(theme_cowplot(font_size = 12))

# BASELINE PARAMETER VALUES =======================================================

p <- 0.15
q <- 0.2
comm <- 550
old.beta1 <- 0.14102773
delta1 <- 0.34
beta2 <- 0.01096798
delta2 <- 0.167 
sigma <- 0.51
gamma1n <- 1/30
gamma1c <- 1/365
gamma2n <- 1/30
gamma2c <- 1/365
gammaD <- 1/5
alphaIL1 <- 0.009
alphaIH1 <- 0.009
alphaIL2 <- 0.009/5
alphaIH2 <- 0.009/5
rSL <- 1/4.13 
rH <- 5.2/4.13 # ~1.259
rI <- 1
rD <- 0.5 
rhoH <- 3.05
rhoI <- 1
rhoD <- 2 
int <- 0 # not in sensitivity analysis
time.int <- 30000 # not in sensitivity analysis 
beta.int <- 0.30
gamma.intL <- 1/28
gamma.intH <- 1/28
decol <- 0 # not in sensitivity analysis 

# need to run the model with no infecteds to get proportions at equilibrium
N1.t0 <- 300
N2.t0 <- N1.t0*comm 
N1c.t0 <- q*N1.t0
N1n.t0 <- (1-q)*N1.t0
N2c.t0 <- q*N2.t0
N2n.t0 <- (1-q)*N2.t0

IL1c.t0 <- 0
IH1c.t0 <- 0
DL1c.t0 <- 0
DH1c.t0 <- 0
IL2c.t0 <- 0
IH2c.t0 <- 0
DL2c.t0 <- 0
DH2c.t0 <- 0
SL1c.t0 <- N1c.t0*(1-p) - IL1c.t0 - DL1c.t0
SH1c.t0 <- N1c.t0*(p) - IH1c.t0 - DH1c.t0
SL2c.t0 <- N2c.t0*(1-p) - IL2c.t0 - DL2c.t0
SH2c.t0 <- N2c.t0*(p) - IH2c.t0 - DH2c.t0

IL1n.t0 <- 0
IH1n.t0 <- 0
DL1n.t0 <- 0
DH1n.t0 <- 0
IL2n.t0 <- 0
IH2n.t0 <- 0
DL2n.t0 <- 0
DH2n.t0 <- 0
SL1n.t0 <- N1n.t0*(1-p) - IL1n.t0 - DL1n.t0
SH1n.t0 <- N1n.t0*(p) - IH1n.t0 - DH1n.t0
SL2n.t0 <- N2n.t0*(1-p) - IL2n.t0 - DL2n.t0
SH2n.t0 <- N2n.t0*(p) - IH2n.t0 - DH2n.t0

Total.t0 <- sum(c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
                  SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0,
                  SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
                  SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0))

# run to find equilibrium with no infecteds
modpars <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
             gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
             rH, rI, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = 0, 
             gamma.intL, gamma.intH, decol = 0)
mod.init <- c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
              SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0, 
              SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
              SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0, Total.t0)
names(mod.init) <- c("SL1c", "SH1c", "IL1c", "IH1c", "DL1c", "DH1c",
                     "SL2c", "SH2c", "IL2c", "IH2c", "DL2c", "DH2c",
                     "SL1n", "SH1n", "IL1n", "IH1n", "DL1n", "DH1n",
                     "SL2n", "SH2n", "IL2n", "IH2n", "DL2n", "DH2n",
                     "total")
uninf.equil <- runsteady(y = mod.init, func = modsim.dyn, parms = modpars, times = c(0, 10000000))

# size of hospital, community, etc
run_N1 <- sum(uninf.equil$y[c(1:6, 13:18)])
run_N2 <- sum(uninf.equil$y[c(7:12, 19:24)])

# size by risk group and location
run_NL1 <- sum(uninf.equil$y[c("SL1c", "IL1c", "DL1c", "SL1n", "IL1n", "DL1n")])
run_NH1 <- sum(uninf.equil$y[c("SH1c", "IH1c", "DH1c", "SH1n", "IH1n", "DH1n")])
run_NL2 <- sum(uninf.equil$y[c("SL2c", "IL2c", "DL2c", "SL2n", "IL2n", "DL2n")])
run_NH2 <- sum(uninf.equil$y[c("SH2c", "IH2c", "DH2c", "SH2n", "IH2n", "DH2n")])

# size by carrier type and location
run_N1c <- sum(uninf.equil$y[c("SL1c", "IL1c", "DL1c", "SH1c", "IH1c", "DH1c")])
run_N2c <- sum(uninf.equil$y[c("SL2c", "IL2c", "DL2c", "SH2c", "IH2c", "DH2c")])

# proportion high risk (p) and proportion recently discharged (q)
p1 <- run_NH1/run_N1
p2 <- run_NH2/run_N2
q1 <- run_N1c/run_N1
q2 <- run_N2c/run_N2


# fit to normal prev ================================================================

NL2.eq <- sum(uninf.equil$y[c(7, 9, 11, 19, 21, 23)])
NH2.eq <- sum(uninf.equil$y[c(8, 10, 12, 20, 22, 24)])
NL1.eq <- sum(uninf.equil$y[c(1, 3, 5, 13, 15, 17)])


N1.t0 <- 300
N2.t0 <- N1.t0*comm 
N1c.t0 <- q1*N1.t0
N1n.t0 <- (1-q1)*N1.t0
N2c.t0 <- q2*N2.t0
N2n.t0 <- (1-q2)*N2.t0

IL1c.t0 <- 1
IH1c.t0 <- 1
DL1c.t0 <- 0
DH1c.t0 <- 0
IL2c.t0 <- 1
IH2c.t0 <- 1
DL2c.t0 <- 0
DH2c.t0 <- 0
SL1c.t0 <- N1c.t0*(1-p1) - IL1c.t0 - DL1c.t0
SH1c.t0 <- N1c.t0*(p1) - IH1c.t0 - DH1c.t0
SL2c.t0 <- N2c.t0*(1-p2) - IL2c.t0 - DL2c.t0
SH2c.t0 <- N2c.t0*(p2) - IH2c.t0 - DH2c.t0

IL1n.t0 <- 1
IH1n.t0 <- 1
DL1n.t0 <- 0
DH1n.t0 <- 0
IL2n.t0 <- 1
IH2n.t0 <- 1
DL2n.t0 <- 0
DH2n.t0 <- 0
SL1n.t0 <- N1n.t0*(1-p1) - IL1n.t0 - DL1n.t0
SH1n.t0 <- N1n.t0*(p1) - IH1n.t0 - DH1n.t0
SL2n.t0 <- N2n.t0*(1-p2) - IL2n.t0 - DL2n.t0
SH2n.t0 <- N2n.t0*(p2) - IH2n.t0 - DH2n.t0

Total.t0 <- sum(c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
                  SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0,
                  SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
                  SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0))


# fit the transmission parameters
params0 <- c(log(0.1), log(0.01))  # initial guess, transformed so range is (0, +inf)
fit2 <- optim(params0, sse.rs.prev) # fit
beta.fit <- exp(fit2$par) # back transform

old.beta1 <- beta.fit[1]
beta2 <- beta.fit[2]

# equilibrium with no intervention
modpars <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
             gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
             rH, rI, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = 0, 
             gamma.intL = 1/28, gamma.intH = 1/28, decol = 0)
mod.init <- c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
              SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0, 
              SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
              SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0, Total.t0)
names(mod.init) <- c("SL1c", "SH1c", "IL1c", "IH1c", "DL1c", "DH1c",
                     "SL2c", "SH2c", "IL2c", "IH2c", "DL2c", "DH2c",
                     "SL1n", "SH1n", "IL1n", "IH1n", "DL1n", "DH1n",
                     "SL2n", "SH2n", "IL2n", "IH2n", "DL2n", "DH2n",
                     "total")
no.int.equil <- runsteady(y = mod.init, func = modsim.dyn, parms = modpars, times = c(0, 10000000))

eqprev.hosp.noint <- sum(no.int.equil$y[c(3, 4, 5, 6, 15, 16, 17, 18)])/N1.t0
eqprev.comm.noint <- sum(no.int.equil$y[c(9, 10, 11, 12, 21, 22, 23, 24)])/N2.t0
eqprev.hosp.noint; eqprev.comm.noint

infecteds <- sum(no.int.equil$y[c(3, 4, 5, 6, 15, 16, 17, 18, 9, 10, 11, 12, 21, 22, 23, 24 )])
infecteds

NL2.eq <- sum(no.int.equil$y[c(7, 9, 11, 19, 21, 23)])
NH2.eq <- sum(no.int.equil$y[c(8, 10, 12, 20, 22, 24)])
NL1.eq <- sum(no.int.equil$y[c(1, 3, 5, 13, 15, 17)])


# equilibrium prevalence by age group
eqprev.comm.low <- sum(no.int.equil$y[c(9, 11, 21, 23)])/NL2.eq
eqprev.comm.high <- sum(no.int.equil$y[c(10, 12, 22, 24)])/NH2.eq

eqprev.comm.high # compare to 3.1% from Gorwitz 2008 JID
eqprev.comm.low # compare to 1.1-1.3% from Gorwitz 2008 JID




# PARAMETER FITTING WITH DIFFRENT ASSUMPTIONS ================================================

# (1) set up vectors for each assumption

homogen <- c(sigma = 1, delta1 = 0, delta2 = 0, rSL = 1/4.5, rSH = 1/4.5, 
            rIL = 1/4.5, rIH = 1/4.5, rhoH = 1, rhoI = 1, gamma1n = 1/278, 
            gamma1c = 1/278, gamma2n = 1/278, gamma2c = 1/278)

heterogen <- c(sigma = 1, delta1 = 0, delta2 = 0, rSL = 1/4.5, rSH = 1/4.5, 
            rIL = 1/4.5, rIH = 1/4.5, rhoH = 1, rhoI = 1, gamma1n = 1/30, 
            gamma1c = 1/365, gamma2n = 1/30, gamma2c = 1/365)

baseline <- c(sigma = 0.51, delta1 = 0.35, delta2 = 0.167, rSL = 1/4.13, rSH = 1/5.2, 
              rIL = 1/4.13, rIH = 1/5.2, rhoH = 3, rhoI = 1, gamma1n = 1/30, 
              gamma1c = 1/365, gamma2n = 1/30, gamma2c = 1/365)  

inf_los <- c(sigma = 1, delta1 = 0, delta2 = 0, rSL = 1/4.466, rSH = 1/4.466, 
             rIL = 1/5.466, rIH = 1/5.466, rhoH = 1, rhoI = 1, gamma1n = 1/30, 
             gamma1c = 1/365, gamma2n = 1/30, gamma2c = 1/365)  

inf_adm <- c(sigma = 1, delta1 = 0, delta2 = 0, rSL = 1/4.5, rSH = 1/4.5, 
             rIL = 1/4.5, rIH = 1/4.5, rhoH = 1, rhoI = 1.3, gamma1n = 1/30, 
             gamma1c = 1/365, gamma2n = 1/30, gamma2c = 1/365) 

inf_both <- c(sigma = 1, delta1 = 0, delta2 = 0, rSL = 1/4.466, rSH = 1/4.466, 
              rIL = 1/5.466, rIH = 1/5.466, rhoH = 1, rhoI = 1.3, gamma1n = 1/30, 
              gamma1c = 1/365, gamma2n = 1/30, gamma2c = 1/365) 

baseline_wt <- c(sigma = 0.51, delta1 = 0.35, delta2 = 0.167, rSL = 1/2.23, rSH = 1/5.2, 
              rIL = 1/2.39, rIH = 1/5.2, rhoH = 3.055, rhoI = 1, gamma1n = 1/30, 
              gamma1c = 1/365, gamma2n = 1/30, gamma2c = 1/365) 

# updated
homogen <- c(sigma = 1, delta1 = 0, delta2 = 0, rSL = 1/4.5, rH = 1, 
             rI = 1, rhoH = 1, rhoI = 1, gamma1n = 1/278, 
             gamma1c = 1/278, gamma2n = 1/278, gamma2c = 1/278)

heterogen <- c(sigma = 1, delta1 = 0, delta2 = 0, rSL = 1/4.5, rH = 1, 
               rI = 1, rhoH = 1, rhoI = 1, gamma1n = 1/30, 
               gamma1c = 1/365, gamma2n = 1/30, gamma2c = 1/365)

baseline <- c(sigma = 0.51, delta1 = 0.35, delta2 = 0.167, rSL = 1/4.13, rH = 5.2/4.13, 
              rI = 1, rhoH = 3, rhoI = 1, gamma1n = 1/30, 
              gamma1c = 1/365, gamma2n = 1/30, gamma2c = 1/365)  

inf_los <- c(sigma = 1, delta1 = 0, delta2 = 0, rSL = 1/4.466, rH = 1, 
             rI = 5.466/4.466, rhoH = 1, rhoI = 1, gamma1n = 1/30, 
             gamma1c = 1/365, gamma2n = 1/30, gamma2c = 1/365)  

inf_adm <- c(sigma = 1, delta1 = 0, delta2 = 0, rSL = 1/4.5, rH = 1, 
             rI = 1, rhoH = 1, rhoI = 1.3, gamma1n = 1/30, 
             gamma1c = 1/365, gamma2n = 1/30, gamma2c = 1/365) 

inf_both <- c(sigma = 1, delta1 = 0, delta2 = 0, rSL = 1/4.466, rH = 1, 
              rI = 5.466/4.466, rhoH = 1, rhoI = 1.3, gamma1n = 1/30, 
              gamma1c = 1/365, gamma2n = 1/30, gamma2c = 1/365) 


assumptions_df <- as.data.frame(rbind(homogen, heterogen, baseline, inf_los, inf_adm, inf_both))
assumptions_df$assumption <- rownames(assumptions_df)
rownames(assumptions_df) <- NULL
assumptions_df[, "p1"] <- NA
assumptions_df[, "p2"] <- NA
assumptions_df[, "q1"] <- NA
assumptions_df[, "q2"] <- NA
assumptions_df[, "old.beta1"] <- NA
assumptions_df[, "beta2"] <- NA
assumptions_df[, "eqprev.hosp.noint"] <- NA
assumptions_df[, "eqprev.comm.noint"] <- NA

previous.fit <- read.csv("assumptions_beta_fiting.csv")
# old.beta1 <- previous.fit[1,19]
# beta2 <- previous.fit[1,20]
# previous.fit[1,19] <- 0.127654703
# previous.fit[1,20] <- 0.003116824

# (2) loop to fit betas for different assumptions

# previous.fit <- read.csv("assumptions_beta_fiting.csv") # get earlier fits

for (i in 1:nrow(assumptions_df)) {
  
  print(i)
  
  sigma <- assumptions_df[i, "sigma"]  
  delta1 <- assumptions_df[i, "delta1"]
  delta2 <- assumptions_df[i, "delta2"]
  rSL <- assumptions_df[i, "rSL"] 
  rH <- assumptions_df[i, "rH"] 
  rI <- assumptions_df[i, "rI"] 
  rhoH <- assumptions_df[i, "rhoH"] 
  rhoI <- assumptions_df[i, "rhoI"] 
  gamma1n <- assumptions_df[i, "gamma1n"]   
  gamma1c <- assumptions_df[i, "gamma1c"]  
  gamma2n <- assumptions_df[i, "gamma2n"]  
  gamma2c <- assumptions_df[i, "gamma2c"]  
  
  
  # need to run the model with no infecteds to get proportions at equilibrium
  N1.t0 <- 300
  N2.t0 <- N1.t0*comm 
  N1c.t0 <- q*N1.t0
  N1n.t0 <- (1-q)*N1.t0
  N2c.t0 <- q*N2.t0
  N2n.t0 <- (1-q)*N2.t0
  
  IL1c.t0 <- 0
  IH1c.t0 <- 0
  DL1c.t0 <- 0
  DH1c.t0 <- 0
  IL2c.t0 <- 0
  IH2c.t0 <- 0
  DL2c.t0 <- 0
  DH2c.t0 <- 0
  SL1c.t0 <- N1c.t0*(1-p) - IL1c.t0 - DL1c.t0
  SH1c.t0 <- N1c.t0*(p) - IH1c.t0 - DH1c.t0
  SL2c.t0 <- N2c.t0*(1-p) - IL2c.t0 - DL2c.t0
  SH2c.t0 <- N2c.t0*(p) - IH2c.t0 - DH2c.t0
  
  IL1n.t0 <- 0
  IH1n.t0 <- 0
  DL1n.t0 <- 0
  DH1n.t0 <- 0
  IL2n.t0 <- 0
  IH2n.t0 <- 0
  DL2n.t0 <- 0
  DH2n.t0 <- 0
  SL1n.t0 <- N1n.t0*(1-p) - IL1n.t0 - DL1n.t0
  SH1n.t0 <- N1n.t0*(p) - IH1n.t0 - DH1n.t0
  SL2n.t0 <- N2n.t0*(1-p) - IL2n.t0 - DL2n.t0
  SH2n.t0 <- N2n.t0*(p) - IH2n.t0 - DH2n.t0
  
  Total.t0 <- sum(c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
                    SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0,
                    SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
                    SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0))
  
  # run to find equilibrium with no infecteds
  modpars <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
               gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
               rH, rI, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = 0, 
               gamma.intL, gamma.intH, decol = 0)
  mod.init <- c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
                SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0, 
                SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
                SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0, Total.t0)
  names(mod.init) <- c("SL1c", "SH1c", "IL1c", "IH1c", "DL1c", "DH1c",
                       "SL2c", "SH2c", "IL2c", "IH2c", "DL2c", "DH2c",
                       "SL1n", "SH1n", "IL1n", "IH1n", "DL1n", "DH1n",
                       "SL2n", "SH2n", "IL2n", "IH2n", "DL2n", "DH2n",
                       "total")
  uninf.equil <- runsteady(y = mod.init, func = modsim.dyn, parms = modpars, times = c(0, 10000000))
  
  # size of hospital, community, etc
  run_N1 <- sum(uninf.equil$y[c(1:6, 13:18)])
  run_N2 <- sum(uninf.equil$y[c(7:12, 19:24)])
  
  # size by risk group and location
  run_NL1 <- sum(uninf.equil$y[c("SL1c", "IL1c", "DL1c", "SL1n", "IL1n", "DL1n")])
  run_NH1 <- sum(uninf.equil$y[c("SH1c", "IH1c", "DH1c", "SH1n", "IH1n", "DH1n")])
  run_NL2 <- sum(uninf.equil$y[c("SL2c", "IL2c", "DL2c", "SL2n", "IL2n", "DL2n")])
  run_NH2 <- sum(uninf.equil$y[c("SH2c", "IH2c", "DH2c", "SH2n", "IH2n", "DH2n")])
  
  # size by carrier type and location
  run_N1c <- sum(uninf.equil$y[c("SL1c", "IL1c", "DL1c", "SH1c", "IH1c", "DH1c")])
  run_N2c <- sum(uninf.equil$y[c("SL2c", "IL2c", "DL2c", "SH2c", "IH2c", "DH2c")])
  
  # proportion high risk (p) and proportion recently discharged (q)
  p1 <- run_NH1/run_N1
  p2 <- run_NH2/run_N2
  q1 <- run_N1c/run_N1
  q2 <- run_N2c/run_N2
  
  N1.t0 <- 300
  N2.t0 <- N1.t0*comm 
  N1c.t0 <- q1*N1.t0
  N1n.t0 <- (1-q1)*N1.t0
  N2c.t0 <- q2*N2.t0
  N2n.t0 <- (1-q2)*N2.t0
  
  IL1c.t0 <- 1
  IH1c.t0 <- 1
  DL1c.t0 <- 0
  DH1c.t0 <- 0
  IL2c.t0 <- 1
  IH2c.t0 <- 1
  DL2c.t0 <- 0
  DH2c.t0 <- 0
  SL1c.t0 <- N1c.t0*(1-p1) - IL1c.t0 - DL1c.t0
  SH1c.t0 <- N1c.t0*(p1) - IL1c.t0 - DL1c.t0
  SL2c.t0 <- N2c.t0*(1-p2) - IL2c.t0 - DL2c.t0
  SH2c.t0 <- N2c.t0*(p2) - IL2c.t0 - DL2c.t0
  
  IL1n.t0 <- 1
  IH1n.t0 <- 1
  DL1n.t0 <- 0
  DH1n.t0 <- 0
  IL2n.t0 <- 1
  IH2n.t0 <- 1
  DL2n.t0 <- 0
  DH2n.t0 <- 0
  SL1n.t0 <- N1n.t0*(1-p1) - IL1n.t0 - DL1n.t0
  SH1n.t0 <- N1n.t0*(p1) - IL1n.t0 - DL1n.t0
  SL2n.t0 <- N2n.t0*(1-p2) - IL2n.t0 - DL2n.t0
  SH2n.t0 <- N2n.t0*(p2) - IL2n.t0 - DL2n.t0
  
  Total.t0 <- sum(c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
                    SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0,
                    SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
                    SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0))
  
  
  # fit transmission parameters
  #params0 <- c(log(0.14), log(0.01)) # use this instead if no previous fits
  params0 <- c(log(previous.fit[i, "old.beta1"]), log(previous.fit[i, "beta2"]))
  fit2 <- optim(params0, sse.rs.prev) # fit
  beta.fit <- exp(fit2$par)
  
  # find equilibrium prevalences with new fit parameters
  old.beta1 <- beta.fit[1]
  beta2 <- beta.fit[2]
  
  # equilibrium with no intervention
  modpars <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
               gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
               rH, rI, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = 0, 
               gamma.intL, gamma.intH, decol = 0)
  mod.init <- c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
                SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0, 
                SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
                SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0, Total.t0)
  names(mod.init) <- c("SL1c", "SH1c", "IL1c", "IH1c", "DL1c", "DH1c",
                       "SL2c", "SH2c", "IL2c", "IH2c", "DL2c", "DH2c",
                       "SL1n", "SH1n", "IL1n", "IH1n", "DL1n", "DH1n",
                       "SL2n", "SH2n", "IL2n", "IH2n", "DL2n", "DH2n",
                       "total")
  no.int.equil <- runsteady(y = mod.init, func = modsim.dyn, parms = modpars, times = c(0, 10000000))
  
  eqprev.hosp.noint <- sum(no.int.equil$y[c(3, 4, 5, 6, 15, 16, 17, 18)])/N1.t0
  eqprev.comm.noint <- sum(no.int.equil$y[c(9, 10, 11, 12, 21, 22, 23, 24)])/N2.t0
  
  
  assumptions_df[i, "p1"] <- p1
  assumptions_df[i, "p2"] <- p2
  assumptions_df[i, "q1"] <- q1
  assumptions_df[i, "q2"] <- q2
  assumptions_df[i, "old.beta1"] <- old.beta1
  assumptions_df[i, "beta2"] <- beta2
  assumptions_df[i, "eqprev.hosp.noint"] <- eqprev.hosp.noint
  assumptions_df[i, "eqprev.comm.noint"] <- eqprev.comm.noint
  
}

# save data
write.csv(assumptions_df, file = "assumptions_beta_fiting2.csv", row.names = FALSE)


# SAVED BETA VALUES FOR DIFFERENT ASSUMPTIONS ASSUMPTIONS ==================================
assumptions_df <- read.csv("assumptions_beta_fiting2.csv")


# EFFECT OF VARYING INTERVENTION EFFICACY--TRANSMISSION =====================================

beta.int.vec <- seq(0, 1, length.out = 25)

beta.data <- as.data.frame(matrix(rep(NA, 8*25*6), nrow = 25*6))
names(beta.data) <- c("beta.int", "R0", "hospital", "community", "cases.averted", 
                                         "baseline.cases", "prop.cases.averted" ,"assumption")
j <- 1

for (i in 1:nrow(assumptions_df)) {
  
  sigma <- assumptions_df[i, "sigma"] 
  delta1 <- assumptions_df[i, "delta1"]
  delta2 <- assumptions_df[i, "delta2"]
  rSL <- assumptions_df[i, "rSL"] 
  rH <- assumptions_df[i, "rH"] 
  rI <- assumptions_df[i, "rI"] 
  rhoH <- assumptions_df[i, "rhoH"] 
  rhoI <- assumptions_df[i, "rhoI"] 
  p1 <- assumptions_df[i, "p1"]
  p2 <- assumptions_df[i, "p2"] 
  q1 <- assumptions_df[i, "q1"]
  q2 <- assumptions_df[i, "q2"]  
  old.beta1 <- assumptions_df[i, "old.beta1"]
  beta2 <- assumptions_df[i, "beta2"]  
  gamma1n <- assumptions_df[i, "gamma1n"]   
  gamma1c <- assumptions_df[i, "gamma1c"]  
  gamma2n <- assumptions_df[i, "gamma2n"]  
  gamma2c <- assumptions_df[i, "gamma2c"] 
  
  for (h in 1:length(beta.int.vec)) {
    
    new.beta.int <- beta.int.vec[h]
    
    # calculate R0 for DFE
    N1.t0 <- 300
    N2.t0 <- N1.t0*comm 
    N1c.t0 <- q1*N1.t0
    N1n.t0 <- (1-q1)*N1.t0
    N2c.t0 <- q2*N2.t0
    N2n.t0 <- (1-q2)*N2.t0
    
    IL1c.t0 <- 0
    IH1c.t0 <- 0
    DL1c.t0 <- 0
    DH1c.t0 <- 0
    IL2c.t0 <- 0
    IH2c.t0 <- 0
    DL2c.t0 <- 0
    DH2c.t0 <- 0
    SL1c.t0 <- N1c.t0*(1-p1) - IL1c.t0 - DL1c.t0
    SH1c.t0 <- N1c.t0*(p1) - IH1c.t0 - DH1c.t0
    SL2c.t0 <- N2c.t0*(1-p2) - IL2c.t0 - DL2c.t0
    SH2c.t0 <- N2c.t0*(p2) - IH2c.t0 - DH2c.t0
    
    IL1n.t0 <- 0
    IH1n.t0 <- 0
    DL1n.t0 <- 0
    DH1n.t0 <- 0
    IL2n.t0 <- 0
    IH2n.t0 <- 0
    DL2n.t0 <- 0
    DH2n.t0 <- 0
    SL1n.t0 <- N1n.t0*(1-p1) - IL1n.t0 - DL1n.t0
    SH1n.t0 <- N1n.t0*(p1) - IH1n.t0 - DH1n.t0
    SL2n.t0 <- N2n.t0*(1-p2) - IL2n.t0 - DL2n.t0
    SH2n.t0 <- N2n.t0*(p2) - IH2n.t0 - DH2n.t0
    
    Total.t0 <- sum(c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
                      SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0,
                      SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
                      SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0))
    
    # set new hospital transmission parameter and calc R0 with it
    new.beta1 <- old.beta1*(1 - new.beta.int)
    R0etc <- calc_R0(comm, q1, q2, p1, p2, rSL, rH, rI, rD, rhoH, rhoI, rhoD, old.beta1 = new.beta1, int, delta1, beta2, sigma, delta2, 
                     gamma1c, alphaIL1, gamma1n, alphaIH1, gamma2c, alphaIL2, gamma2n, alphaIH2, gammaD)
    
    # add infectious indivs for simulations
    N1.t0 <- 300
    N2.t0 <- N1.t0*comm 
    N1c.t0 <- q1*N1.t0
    N1n.t0 <- (1-q1)*N1.t0
    N2c.t0 <- q2*N2.t0
    N2n.t0 <- (1-q2)*N2.t0
    
    IL1c.t0 <- 1
    IH1c.t0 <- 1
    DL1c.t0 <- 0
    DH1c.t0 <- 0
    IL2c.t0 <- 1
    IH2c.t0 <- 1
    DL2c.t0 <- 0
    DH2c.t0 <- 0
    SL1c.t0 <- N1c.t0*(1-p1) - IL1c.t0 - DL1c.t0
    SH1c.t0 <- N1c.t0*(p1) - IH1c.t0 - DH1c.t0
    SL2c.t0 <- N2c.t0*(1-p2) - IL2c.t0 - DL2c.t0
    SH2c.t0 <- N2c.t0*(p2) - IH2c.t0 - DH2c.t0
    
    IL1n.t0 <- 1
    IH1n.t0 <- 1
    DL1n.t0 <- 0
    DH1n.t0 <- 0
    IL2n.t0 <- 1
    IH2n.t0 <- 1
    DL2n.t0 <- 0
    DH2n.t0 <- 0
    SL1n.t0 <- N1n.t0*(1-p1) - IL1n.t0 - DL1n.t0
    SH1n.t0 <- N1n.t0*(p1) - IH1n.t0 - DH1n.t0
    SL2n.t0 <- N2n.t0*(1-p2) - IL2n.t0 - DL2n.t0
    SH2n.t0 <- N2n.t0*(p2) - IH2n.t0 - DH2n.t0
    
    Total.t0 <- sum(c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
                      SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0,
                      SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
                      SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0))
    
    # 1 No intervention--runsteady to find equilibrium
    # equilibrium with no intervention
    modpars <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
                 gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
                 rH, rI, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = 0, 
                 gamma.intL, gamma.intH, decol = 0)
    mod.init <- c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
                  SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0, 
                  SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
                  SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0, Total.t0)
    names(mod.init) <- c("SL1c", "SH1c", "IL1c", "IH1c", "DL1c", "DH1c",
                         "SL2c", "SH2c", "IL2c", "IH2c", "DL2c", "DH2c",
                         "SL1n", "SH1n", "IL1n", "IH1n", "DL1n", "DH1n",
                         "SL2n", "SH2n", "IL2n", "IH2n", "DL2n", "DH2n",
                         "total")
    no.int.equil <- runsteady(y = mod.init, func = modsim.dyn, parms = modpars, times = c(0, 10000000))
    equil.init <- unname(no.int.equil$y) # save like this so it's easier to use
    
    eqprev.hosp.noint <- sum(no.int.equil$y[c(3, 4, 5, 6, 15, 16, 17, 18)])/N1.t0
    eqprev.comm.noint <- sum(no.int.equil$y[c(9, 10, 11, 12, 21, 22, 23, 24)])/N2.t0
    
    
    # 2 No intervention--deSolve to find 5 year incidence
    modpars <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
                 gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
                 rH, rI, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = 0, 
                 gamma.intL, gamma.intH, decol = 0)
    mod.init <- c(SL1c.t0 = equil.init[1], SH1c.t0 = equil.init[2], IL1c.t0 = equil.init[3], 
                  IH1c.t0 = equil.init[4], DL1c.t0 = equil.init[5], DH1c.t0 = equil.init[6], 
                  SL2c.t0 = equil.init[7], SH2c.t0 = equil.init[8], IL2c.t0 = equil.init[9], 
                  IH2c.t0 = equil.init[10], DL2c.t0 = equil.init[11], DH2c.t0 = equil.init[12], 
                  SL1n.t0 = equil.init[13], SH1n.t0 = equil.init[14], IL1n.t0 = equil.init[15], 
                  IH1n.t0 = equil.init[16], DL1n.t0 = equil.init[17], DH1n.t0 = equil.init[18], 
                  SL2n.t0 = equil.init[19], SH2n.t0 = equil.init[20], IL2n.t0 = equil.init[21], 
                  IH2n.t0 = equil.init[22], DL2n.t0 = equil.init[23], DH2n.t0 = equil.init[24],
                  Total.t0 = equil.init[25])
    mod.t <- seq(0, 1825, by=1) 
    mod.sol.noint <- as.data.frame(lsoda(mod.init, mod.t, modsim.dyn, modpars))
    
    baseline.trans <- sum(mod.sol.noint[, 27:34]) # incidence of transmission (total)
    baseline.trans.hosp <- sum(mod.sol.noint[, c(27, 28, 31, 32)]) # incidence of transmission (hospital)
    baseline.trans.comm <- sum(mod.sol.noint[, c(29, 30, 33, 34)]) # incidence of transmission (community)
    
    baseline.cases <- sum(mod.sol.noint[, 35:42]) # incidence of disease (total)
    baseline.cases.hosp <- sum(mod.sol.noint[, c(35, 36, 39, 40)]) # incidence of disease (hospital)
    baseline.cases.comm <- sum(mod.sol.noint[, c(37, 38, 41, 42)]) # incidence of disease (community)
    
    # 3 Intervention--deSolve to find 5 year incidence
    modpars.int <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
                     gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
                     rH, rI, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = new.beta.int, 
                     gamma.intL, gamma.intH, decol = 0)
    mod.init <- c(SL1c.t0 = equil.init[1], SH1c.t0 = equil.init[2], IL1c.t0 = equil.init[3], 
                  IH1c.t0 = equil.init[4], DL1c.t0 = equil.init[5], DH1c.t0 = equil.init[6], 
                  SL2c.t0 = equil.init[7], SH2c.t0 = equil.init[8], IL2c.t0 = equil.init[9], 
                  IH2c.t0 = equil.init[10], DL2c.t0 = equil.init[11], DH2c.t0 = equil.init[12], 
                  SL1n.t0 = equil.init[13], SH1n.t0 = equil.init[14], IL1n.t0 = equil.init[15], 
                  IH1n.t0 = equil.init[16], DL1n.t0 = equil.init[17], DH1n.t0 = equil.init[18], 
                  SL2n.t0 = equil.init[19], SH2n.t0 = equil.init[20], IL2n.t0 = equil.init[21], 
                  IH2n.t0 = equil.init[22], DL2n.t0 = equil.init[23], DH2n.t0 = equil.init[24],
                  Total.t0 = equil.init[25])
    mod.t <- seq(0, 1825, by=1) 
    mod.sol.int <- as.data.frame(lsoda(mod.init, mod.t, modsim.dyn, modpars.int))
    
    intervention.trans <- sum(mod.sol.int[, 27:34]) # incidence of transmission (total)
    intervention.trans.hosp <- sum(mod.sol.int[, c(27, 28, 31, 32)]) # incidence of transmission (hospital)
    intervention.trans.comm <- sum(mod.sol.int[, c(29, 30, 33, 34)]) # incidence of transmission (community)
    
    intervention.cases <- sum(mod.sol.int[, 35:42]) # incidence of disease (total)
    intervention.cases.hosp <- sum(mod.sol.int[, c(35, 36, 39, 40)]) # incidence of disease (hospital)
    intervention.cases.comm <- sum(mod.sol.int[, c(37, 38, 41, 42)]) # incidence of disease (community)
    
    cases.averted <- baseline.cases - intervention.cases
    trans.averted <- baseline.trans  - intervention.trans
    prop.cases.averted <- cases.averted/baseline.cases
    
    # 4 Intervention--runsteady to find equilibrium 
    modpars.int <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
                     gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
                     rH, rI, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = new.beta.int, 
                     gamma.intL, gamma.intH, decol = 0)
    mod.init <- c(SL1c.t0 = equil.init[1], SH1c.t0 = equil.init[2], IL1c.t0 = equil.init[3], 
                  IH1c.t0 = equil.init[4], DL1c.t0 = equil.init[5], DH1c.t0 = equil.init[6], 
                  SL2c.t0 = equil.init[7], SH2c.t0 = equil.init[8], IL2c.t0 = equil.init[9], 
                  IH2c.t0 = equil.init[10], DL2c.t0 = equil.init[11], DH2c.t0 = equil.init[12], 
                  SL1n.t0 = equil.init[13], SH1n.t0 = equil.init[14], IL1n.t0 = equil.init[15], 
                  IH1n.t0 = equil.init[16], DL1n.t0 = equil.init[17], DH1n.t0 = equil.init[18], 
                  SL2n.t0 = equil.init[19], SH2n.t0 = equil.init[20], IL2n.t0 = equil.init[21], 
                  IH2n.t0 = equil.init[22], DL2n.t0 = equil.init[23], DH2n.t0 = equil.init[24],
                  Total.t0 = equil.init[25])
    int.equil <- runsteady(y = mod.init, func = modsim.dyn, parms = modpars.int, times = c(0, 10000000))
    
    eqprev.hosp.int <- sum(int.equil$y[c(3, 4, 5, 6, 15, 16, 17, 18)])/N1.t0
    eqprev.comm.int <- sum(int.equil$y[c(9, 10, 11, 12, 21, 22, 23, 24)])/N2.t0
    
    beta.data[j, 1] <- beta.int.vec[h]
    beta.data[j, 2] <- R0etc$R0
    beta.data[j, 3] <- eqprev.hosp.int
    beta.data[j, 4] <- eqprev.comm.int
    beta.data[j, 5] <- cases.averted
    beta.data[j, 6] <- baseline.cases
    beta.data[j, 7] <- prop.cases.averted
    beta.data[j, 8] <- paste(assumptions_df[i, "assumption"])
    
    j <- j + 1
    print(j)
  }
}

# save data
write.csv(beta.data, file = "varying_beta_int.csv", row.names = FALSE)

beta.data$assumption <- as.factor(beta.data$assumption)

# color scheme for plotting
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

# ignore heterogeneous example
main.beta <- beta.data %>% filter(assumption != "heterogen")

# plot raw number of cases averted
beta.int.cases.all <- ggplot(data = main.beta, aes(x = beta.int*100, y = cases.averted, color = assumption)) +
  geom_line(size = 1) +
  ylab("Cases averted") +
  xlab(expression(atop("Transmission intervention", paste("effectiveness (%, ", theta,")")))) +
  ggtitle("") +
  #scale_y_continuous(limits = c(0, 3500), breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500)) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100")) +
  scale_color_manual(breaks = c("homogen", "baseline",  "inf_adm", "inf_los", "inf_both"), 
                     values = cbbPalette,
                     labels = c("Homogeneous", "Age group", 
                                 "Colonized Admission", "Colonized LOS", "Colonized LOS +\nAdmission"),
                     guide_legend(title = "Assumption")) 
beta.int.cases.all 

# plot proportional number of cases averted
beta.int.propcases <- ggplot(data = main.beta, aes(x = beta.int*100, y = prop.cases.averted*100, color = assumption)) +
  geom_line(size = 1) +
  ylab("Percentage of cases averted") +
  xlab(expression(atop("Transmission intervention", paste("effectiveness (%, ", theta,")")))) +
  ggtitle("") +
  scale_y_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +  
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100")) +
  scale_color_manual(breaks = c("homogen", "baseline",  "inf_adm", "inf_los", "inf_both"), 
                     values = cbbPalette,
                     labels = c("Homogeneous", "Risk group", 
                                "Colonized Admission", "Colonized LOS", "Colonized LOS +\nAdmission"),
                     guide_legend(title = "Assumption")) 
beta.int.propcases 

# plot pseudo R0 value
beta.int.R0 <- ggplot(data = main.beta, aes(x = beta.int, y = R0, color = assumption)) +
  geom_line(size = 1) +
  ylab("R0") +
  xlab("Transmission intervention efficacy") +
  ggtitle("pseudo-R0 value with transmission intervention") +
  geom_hline(yintercept = 1, alpha = 0.4, linetype = "dashed") + 
  #scale_y_continuous(limits = c(0, 1500), breaks = c(0, 500, 1000, 1500)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.50", "0.75", "1"), expand = c(0, 0)) +
  theme(legend.position = "right") +
  scale_color_manual(breaks = c("homogen", "baseline",  "inf_adm", "inf_los", "inf_both"), 
                     values = cbbPalette,
                     labels = c("Homogeneous carriage", "Risk group", 
                                "Colonized admission", "Colonized LOS", "Both")) 
beta.int.R0 

# organize data for plotting equil. prevalence by location
tidy.beta.data <- beta.data %>% 
  gather(Location, Prevalence, hospital:community)
tidy.beta.data$assumption <- as.factor(tidy.beta.data$assumption)
tidy.beta.data$Location <- as.factor(tidy.beta.data$Location)

# plot equil. prevalence
beta.int.prev <- ggplot(data = tidy.beta.data, aes(x = beta.int, y = Prevalence*100, color = Location)) +
  #geom_hline(yintercept = 0, alpha = 0.4, linetype = "dashed") + 
  geom_line(size = 1) +
  ylab("prevalence (%)") +
  xlab("Transmission intervention efficacy") +
  ggtitle("After an intervention, what is the equilibrium prevalence in each location?") +
  scale_color_manual(values = c("#1b9e77", "#7570b3"),
                     labels = c("Community", "Hospital")) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.50", "0.75", "1")) +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) + 
  facet_wrap(~assumption, ncol = 6) +
  theme(legend.position = "bottom")
beta.int.prev 

# save plots as pdfs
# pdf(file = "Cases_by_assumption_betavec.pdf", width = 5.5, height = 4, useDingbats=F)
# beta.int.cases.all
# dev.off()
# 
# pdf(file = "R0_cases_by_assumption_betavec.pdf", width = 5, height = 4, useDingbats=F)
# beta.int.cases.all
# beta.int.R0
# dev.off()
# 
# pdf(file = "eqprev_assumption_betavec.pdf", width = 8, height = 4, useDingbats=F)
# beta.int.prev
# dev.off()


# EFFECT OF VARYING INTERVENTION EFFICACY--DECOLONIZATION =====================================

decol.int.vec <- seq(1/5, 1/365, length.out = 25) # baseline carriage length = 275 days

decol.data <- as.data.frame(matrix(rep(NA, 8*175), nrow = 175))
names(decol.data) <- names(decol.data) <-c("decol.rate",  "R0", "hospital", "community", "cases.averted", 
                                           "baseline.cases", "prop.cases.averted" ,"assumption")

j <- 1

for (i in 1:nrow(assumptions_df)) {
  
  sigma <- assumptions_df[i, "sigma"] 
  delta1 <- assumptions_df[i, "delta1"]
  delta2 <- assumptions_df[i, "delta2"]
  rSL <- assumptions_df[i, "rSL"] 
  rH <- assumptions_df[i, "rH"] 
  rI <- assumptions_df[i, "rI"] 
  rhoH <- assumptions_df[i, "rhoH"] 
  rhoI <- assumptions_df[i, "rhoI"] 
  p1 <- assumptions_df[i, "p1"]
  p2 <- assumptions_df[i, "p2"] 
  q1 <- assumptions_df[i, "q1"]
  q2 <- assumptions_df[i, "q2"]  
  old.beta1 <- assumptions_df[i, "old.beta1"]
  beta2 <- assumptions_df[i, "beta2"]  
  gamma1n <- assumptions_df[i, "gamma1n"]   
  gamma1c <- assumptions_df[i, "gamma1c"]  
  gamma2n <- assumptions_df[i, "gamma2n"]  
  gamma2c <- assumptions_df[i, "gamma2c"] 
  
  for (h in 1:length(decol.int.vec)) {
    decol.int <- decol.int.vec[h] # vary intervention efficacy
    
    gamma.intH <- decol.int
    gamma.intL <- decol.int
    
    # calculate R0 for DFE
    N1.t0 <- 300
    N2.t0 <- N1.t0*comm 
    N1c.t0 <- q1*N1.t0
    N1n.t0 <- (1-q1)*N1.t0
    N2c.t0 <- q2*N2.t0
    N2n.t0 <- (1-q2)*N2.t0
    
    IL1c.t0 <- 0
    IH1c.t0 <- 0
    DL1c.t0 <- 0
    DH1c.t0 <- 0
    IL2c.t0 <- 0
    IH2c.t0 <- 0
    DL2c.t0 <- 0
    DH2c.t0 <- 0
    SL1c.t0 <- N1c.t0*(1-p1) - IL1c.t0 - DL1c.t0
    SH1c.t0 <- N1c.t0*(p1) - IH1c.t0 - DH1c.t0
    SL2c.t0 <- N2c.t0*(1-p2) - IL2c.t0 - DL2c.t0
    SH2c.t0 <- N2c.t0*(p2) - IH2c.t0 - DH2c.t0
    
    IL1n.t0 <- 0
    IH1n.t0 <- 0
    DL1n.t0 <- 0
    DH1n.t0 <- 0
    IL2n.t0 <- 0
    IH2n.t0 <- 0
    DL2n.t0 <- 0
    DH2n.t0 <- 0
    SL1n.t0 <- N1n.t0*(1-p1) - IL1n.t0 - DL1n.t0
    SH1n.t0 <- N1n.t0*(p1) - IH1n.t0 - DH1n.t0
    SL2n.t0 <- N2n.t0*(1-p2) - IL2n.t0 - DL2n.t0
    SH2n.t0 <- N2n.t0*(p2) - IH2n.t0 - DH2n.t0
    
    Total.t0 <- sum(c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
                      SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0,
                      SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
                      SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0))
    
    # pseudo-R0 doesn't work with decolonization anymore
    # new.beta1 <- old.beta1*(1 - new.beta.int)
    # R0etc <- calc_R0(comm, q1, q2, p1, p2, rD, rIL, rIH, rSL, rSH, rhoH, rhoI, rhoD, old.beta1, int, delta1, beta2, sigma, delta2, 
    #                  gamma1c = gamma1c + gamma.intH, alphaIL1, gamma1n = gamma.intH, alphaIH1, gamma2c, alphaIL2, gamma2n, alphaIH2, gammaD)
    
    # add infectious indivs for simulations
    N1.t0 <- 300
    N2.t0 <- N1.t0*comm 
    N1c.t0 <- q1*N1.t0
    N1n.t0 <- (1-q1)*N1.t0
    N2c.t0 <- q2*N2.t0
    N2n.t0 <- (1-q2)*N2.t0
    
    IL1c.t0 <- 1
    IH1c.t0 <- 1
    DL1c.t0 <- 0
    DH1c.t0 <- 0
    IL2c.t0 <- 1
    IH2c.t0 <- 1
    DL2c.t0 <- 0
    DH2c.t0 <- 0
    SL1c.t0 <- N1c.t0*(1-p1) - IL1c.t0 - DL1c.t0
    SH1c.t0 <- N1c.t0*(p1) - IH1c.t0 - DH1c.t0
    SL2c.t0 <- N2c.t0*(1-p2) - IL2c.t0 - DL2c.t0
    SH2c.t0 <- N2c.t0*(p2) - IH2c.t0 - DH2c.t0
    
    IL1n.t0 <- 1
    IH1n.t0 <- 1
    DL1n.t0 <- 0
    DH1n.t0 <- 0
    IL2n.t0 <- 1
    IH2n.t0 <- 1
    DL2n.t0 <- 0
    DH2n.t0 <- 0
    SL1n.t0 <- N1n.t0*(1-p1) - IL1n.t0 - DL1n.t0
    SH1n.t0 <- N1n.t0*(p1) - IH1n.t0 - DH1n.t0
    SL2n.t0 <- N2n.t0*(1-p2) - IL2n.t0 - DL2n.t0
    SH2n.t0 <- N2n.t0*(p2) - IH2n.t0 - DH2n.t0
    
    Total.t0 <- sum(c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
                      SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0,
                      SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
                      SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0))
    
    # 1 No intervention--runsteady to find equilibrium
    # equilibrium with no intervention
    modpars <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
                 gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
                 rH, rI, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = 0, 
                 gamma.intL, gamma.intH, decol = 0)
    mod.init <- c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
                  SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0, 
                  SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
                  SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0, Total.t0)
    names(mod.init) <- c("SL1c", "SH1c", "IL1c", "IH1c", "DL1c", "DH1c",
                         "SL2c", "SH2c", "IL2c", "IH2c", "DL2c", "DH2c",
                         "SL1n", "SH1n", "IL1n", "IH1n", "DL1n", "DH1n",
                         "SL2n", "SH2n", "IL2n", "IH2n", "DL2n", "DH2n",
                         "total")
    no.int.equil <- runsteady(y = mod.init, func = modsim.dyn, parms = modpars, times = c(0, 10000000))
    equil.init <- unname(no.int.equil$y) # save like this so it's easier to use
    
    eqprev.hosp.noint <- sum(no.int.equil$y[c(3, 4, 5, 6, 15, 16, 17, 18)])/N1.t0
    eqprev.comm.noint <- sum(no.int.equil$y[c(9, 10, 11, 12, 21, 22, 23, 24)])/N2.t0
    
    # 2 No intervention--deSolve to find 5 year incidence
    modpars <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
                 gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
                 rH, rI, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = 0, 
                 gamma.intL, gamma.intH, decol = 0)
    mod.init <- c(SL1c.t0 = equil.init[1], SH1c.t0 = equil.init[2], IL1c.t0 = equil.init[3], 
                  IH1c.t0 = equil.init[4], DL1c.t0 = equil.init[5], DH1c.t0 = equil.init[6], 
                  SL2c.t0 = equil.init[7], SH2c.t0 = equil.init[8], IL2c.t0 = equil.init[9], 
                  IH2c.t0 = equil.init[10], DL2c.t0 = equil.init[11], DH2c.t0 = equil.init[12], 
                  SL1n.t0 = equil.init[13], SH1n.t0 = equil.init[14], IL1n.t0 = equil.init[15], 
                  IH1n.t0 = equil.init[16], DL1n.t0 = equil.init[17], DH1n.t0 = equil.init[18], 
                  SL2n.t0 = equil.init[19], SH2n.t0 = equil.init[20], IL2n.t0 = equil.init[21], 
                  IH2n.t0 = equil.init[22], DL2n.t0 = equil.init[23], DH2n.t0 = equil.init[24],
                  Total.t0 = equil.init[25])
    mod.t <- seq(0, 1825, by=1) 
    mod.sol.noint <- as.data.frame(lsoda(mod.init, mod.t, modsim.dyn, modpars))
    
    baseline.trans <- sum(mod.sol.noint[, 27:34]) # incidence of transmission (total)
    baseline.trans.hosp <- sum(mod.sol.noint[, c(27, 28, 31, 32)]) # incidence of transmission (hospital)
    baseline.trans.comm <- sum(mod.sol.noint[, c(29, 30, 33, 34)]) # incidence of transmission (community)
    
    baseline.cases <- sum(mod.sol.noint[, 35:42]) # incidence of disease (total)
    baseline.cases.hosp <- sum(mod.sol.noint[, c(35, 36, 39, 40)]) # incidence of disease (hospital)
    baseline.cases.comm <- sum(mod.sol.noint[, c(37, 38, 41, 42)]) # incidence of disease (community)
    
    # 3 Intervention--deSolve to find 5 year incidence
    modpars.int <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
                     gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
                     rH, rI, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = 0, 
                     gamma.intL, gamma.intH, decol = 1)
    mod.init <- c(SL1c.t0 = equil.init[1], SH1c.t0 = equil.init[2], IL1c.t0 = equil.init[3], 
                  IH1c.t0 = equil.init[4], DL1c.t0 = equil.init[5], DH1c.t0 = equil.init[6], 
                  SL2c.t0 = equil.init[7], SH2c.t0 = equil.init[8], IL2c.t0 = equil.init[9], 
                  IH2c.t0 = equil.init[10], DL2c.t0 = equil.init[11], DH2c.t0 = equil.init[12], 
                  SL1n.t0 = equil.init[13], SH1n.t0 = equil.init[14], IL1n.t0 = equil.init[15], 
                  IH1n.t0 = equil.init[16], DL1n.t0 = equil.init[17], DH1n.t0 = equil.init[18], 
                  SL2n.t0 = equil.init[19], SH2n.t0 = equil.init[20], IL2n.t0 = equil.init[21], 
                  IH2n.t0 = equil.init[22], DL2n.t0 = equil.init[23], DH2n.t0 = equil.init[24],
                  Total.t0 = equil.init[25])
    mod.t <- seq(0, 1825, by=1) 
    mod.sol.int <- as.data.frame(lsoda(mod.init, mod.t, modsim.dyn, modpars.int))
    
    intervention.trans <- sum(mod.sol.int[, 27:34]) # incidence of transmission (total)
    intervention.trans.hosp <- sum(mod.sol.int[, c(27, 28, 31, 32)]) # incidence of transmission (hospital)
    intervention.trans.comm <- sum(mod.sol.int[, c(29, 30, 33, 34)]) # incidence of transmission (community)
    
    intervention.cases <- sum(mod.sol.int[, 35:42]) # incidence of disease (total)
    intervention.cases.hosp <- sum(mod.sol.int[, c(35, 36, 39, 40)]) # incidence of disease (hospital)
    intervention.cases.comm <- sum(mod.sol.int[, c(37, 38, 41, 42)]) # incidence of disease (community)
    
    cases.averted <- baseline.cases - intervention.cases
    trans.averted <- baseline.trans  - intervention.trans
    prop.cases.averted <- cases.averted/baseline.cases
    
    # 4 Intervention--runsteady to find equilibrium 
    modpars.int <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
                     gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
                     rH, rI, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = 0, 
                     gamma.intL, gamma.intH, decol = 1)
    mod.init <- c(SL1c.t0 = equil.init[1], SH1c.t0 = equil.init[2], IL1c.t0 = equil.init[3], 
                  IH1c.t0 = equil.init[4], DL1c.t0 = equil.init[5], DH1c.t0 = equil.init[6], 
                  SL2c.t0 = equil.init[7], SH2c.t0 = equil.init[8], IL2c.t0 = equil.init[9], 
                  IH2c.t0 = equil.init[10], DL2c.t0 = equil.init[11], DH2c.t0 = equil.init[12], 
                  SL1n.t0 = equil.init[13], SH1n.t0 = equil.init[14], IL1n.t0 = equil.init[15], 
                  IH1n.t0 = equil.init[16], DL1n.t0 = equil.init[17], DH1n.t0 = equil.init[18], 
                  SL2n.t0 = equil.init[19], SH2n.t0 = equil.init[20], IL2n.t0 = equil.init[21], 
                  IH2n.t0 = equil.init[22], DL2n.t0 = equil.init[23], DH2n.t0 = equil.init[24],
                  Total.t0 = equil.init[25])
    int.equil <- runsteady(y = mod.init, func = modsim.dyn, parms = modpars.int, times = c(0, 10000000))
    
    eqprev.hosp.int <- sum(int.equil$y[c(3, 4, 5, 6, 15, 16, 17, 18)])/N1.t0
    eqprev.comm.int <- sum(int.equil$y[c(9, 10, 11, 12, 21, 22, 23, 24)])/N2.t0
    
    decol.data[j, 1] <- decol.int.vec[h]
    decol.data[j, 2] <- NA # R0 calc doesn't work yet with new clearance rate
    decol.data[j, 3] <- eqprev.hosp.int
    decol.data[j, 4] <- eqprev.comm.int
    decol.data[j, 5] <- cases.averted
    decol.data[j, 6] <- baseline.cases
    decol.data[j, 7] <- prop.cases.averted
    decol.data[j, 8] <- paste(assumptions_df[i, "assumption"])
    
    j <- j + 1
    print(j)
  }
}

# save data
write.csv(decol.data, file = "varying_decol_int.csv", row.names = FALSE)

decol.data$assumption <- as.factor(decol.data$assumption)
main.decol <- decol.data %>% filter(assumption != "heterogen" &
                                    assumption != "baseline_wt")

# plot raw number of cases averted
decol.int.cases.all <- ggplot(data = main.decol, aes(x = decol.rate, y = cases.averted, color = assumption)) +
  geom_line(size = 1) +
  ylab("Cases averted") +
  xlab(expression(atop("Decolonization rate", paste("(", gamma[int],")")))) +
  ggtitle("") +
  scale_y_continuous(limits = c(0, 3500), breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500)) +
  scale_x_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                      labels = c("0", "0.05", "0.1", "0.15", "0.20")) +
  scale_color_manual(breaks = c("homogen", "baseline",  "inf_adm", "inf_los", "inf_both"), 
                     values = cbbPalette,
                     labels = c("Homogeneous", "Age group", 
                                "Colonized Admission", "Colonized LOS", "Colonized LOS + Admission"),
                     guide_legend(title = "Assumption")) 
decol.int.cases.all 

# plot proportional amount of cases averted
decol.int.propcases <- ggplot(data = main.decol, aes(x = decol.rate, y = prop.cases.averted*100, color = assumption)) +
  geom_line(size = 1) +
  ylab("Percentage of cases averted") +
  xlab(expression(atop("Decolonization rate", paste("(", gamma[int],")")))) +
  ggtitle("") +
  scale_y_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +
  scale_x_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     labels = c("0", "0.05", "0.1", "0.15", "0.20")) +
  scale_color_manual(breaks = c("homogen", "baseline",  "inf_adm", "inf_los", "inf_both"), 
                     values = cbbPalette,
                     labels = c("Homogeneous", "Age group", 
                                "Colonized Admission", "Colonized LOS", "Colonized LOS + Admission"),
                     guide_legend(title = "Assumption")) 
decol.int.propcases 


# COMBINED TRANSMISSION & DECOLONIZATION INTERVENTION PLOTS ============================

# loading saved data--so we don't have to rerun to change plotting
decol.data <- read.csv("varying_decol_int.csv")
decol.data$assumption <- as.factor(decol.data$assumption)
main.decol <- decol.data %>% filter(assumption != "heterogen" &
                                    assumption != "baseline_wt")

beta.data<- read.csv("varying_beta_int.csv")
beta.data$assumption <- as.factor(beta.data$assumption)
main.beta <- beta.data %>% filter(assumption != "heterogen" &
                                    assumption != "baseline_wt")

# averages
mean.decol <- decol.data %>% 
  group_by(assumption) %>% 
  summarise(mean(cases.averted))

med.decol <- decol.data %>% 
  group_by(assumption) %>% 
  summarise(median(cases.averted))

# color scheme for plotting
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

# plot raw number of cases averted
decol.int.cases.all <- ggplot(data = main.decol, aes(x = decol.rate, y = cases.averted, color = assumption)) +
  geom_line(size = 1) +
  ylab("Cases averted") +
  xlab(expression(atop("Decolonization rate", paste("(1/days, ", gamma[int],")")))) +
  ggtitle("") +
  scale_y_continuous(limits = c(0, 3500), breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500)) +
  scale_x_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     labels = c("0", "0.05", "0.1", "0.15", "0.20")) +
  scale_color_manual(breaks = c("homogen", "baseline",  "inf_adm", "inf_los", "inf_both"), 
                     values = cbbPalette,
                     labels = c("Homogeneous", "Age group", 
                                "Colonized Admission", "Colonized LOS", "Colonized Admission +\nLOS"),
                     guide_legend(title = "Model")) 
decol.int.cases.all 

# plot proportional amount of cases averted
decol.int.propcases <- ggplot(data = main.decol, aes(x = decol.rate, y = prop.cases.averted*100,
                                                     group = assumption, color = assumption)) +
  geom_line(aes(linetype = assumption), size = 0.75) +
  ylab("Percentage of cases averted") +
  xlab(expression(atop("Decolonization rate", paste("(1/days, ", gamma[int],")")))) +
  ggtitle("") +
  scale_y_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +
  scale_x_continuous(limits = c(0, 0.2), breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     labels = c("0", "0.05", "0.1", "0.15", "0.20")) +
  scale_linetype_manual(breaks = c("homogen", "baseline",  "inf_adm", "inf_los", "inf_both"),
                        values = c(1, 1, 2, 3, 4),
                        labels = c("Homogeneous", "Age group", 
                                   "Colonized Admission", "Colonized LOS", "Colonized Admission +\nLOS"),
                        guide_legend(title = "Model")) +
  scale_color_manual(breaks = c("homogen", "baseline",  "inf_adm", "inf_los", "inf_both"), 
                     values = cbbPalette,
                     labels = c("Homogeneous", "Age group", 
                                "Colonized Admission", "Colonized LOS", "Colonized Admission +\nLOS"),
                     guide_legend(title = "Model")) 
decol.int.propcases 

# plot raw number of cases averted
beta.int.cases.all <- ggplot(data = main.beta, aes(x = beta.int*100, y = cases.averted, color = assumption)) +
  geom_line(size = 1) +
  ylab("Cases averted") +
  xlab(expression(atop("Transmission intervention", paste("effectiveness (%, ", theta,")")))) +
  ggtitle("") +
  #scale_y_continuous(limits = c(0, 3500), breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500)) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100")) +
  scale_color_manual(breaks = c("homogen", "baseline",  "inf_adm", "inf_los", "inf_both"), 
                     values = cbbPalette,
                     labels = c("Homogeneous", "Age group", 
                                "Colonized Admission", "Colonized LOS", "Colonized Admission +\nLOS"),
                     guide_legend(title = "Model")) 
beta.int.cases.all 

# plot proportional number of cases averted
beta.int.propcases <- ggplot(data = main.beta, aes(x = beta.int*100, y = prop.cases.averted*100, 
                                                   group = assumption, color = assumption)) +
  geom_line(aes(linetype = assumption), size = 0.75) +
  ylab("Percentage of cases averted") +
  xlab(expression(atop("Transmission intervention", paste("effectiveness (%, ", theta,")")))) +
  ggtitle("") +
  scale_y_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +  
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100")) +
  scale_linetype_manual(breaks = c("homogen", "baseline",  "inf_adm", "inf_los", "inf_both"),
                        values = c(1, 1, 2, 3, 4),
                        labels = c("Homogeneous", "Age group", 
                                   "Colonized Admission", "Colonized LOS", "Colonized Admission +\nLOS"),
                        guide_legend(title = "Model")) +
  scale_color_manual(breaks = c("homogen", "baseline",  "inf_adm", "inf_los", "inf_both"), 
                     values = cbbPalette,
                     labels = c("Homogeneous", "Age group", 
                                "Colonized Admission", "Colonized LOS", "Colonized Admission +\nLOS"),
                     guide_legend(title = "Model")) 
beta.int.propcases 

# raw number of cases plots
plot.row <- plot_grid(beta.int.cases.all + theme(legend.position = "none"), 
                      decol.int.cases.all + theme(legend.position = "none"),
                      labels=c("A", "B"), ncol = 2)
legend <- get_legend(beta.int.cases.all)

plot.full <- plot_grid(plot.row, legend, rel_widths = c(2, .7))

# pdf(file = "Cases_by_assumption.pdf", width = 10, height = 4, useDingbats=F)
# plot.full
# dev.off()

# combine proportional number of cases averted plots
# this is Fig 2 for publication
plot.row.prop <- plot_grid(beta.int.propcases + theme_cowplot(font_size = 12) + theme(legend.position = "none") , 
                      decol.int.propcases + theme_cowplot(font_size = 12) + theme(legend.position = "none") ,
                      labels=c("A", "B"), ncol = 2)
legend.prop <- get_legend(beta.int.propcases + theme_cowplot(font_size = 12))

plot.fullprop <- plot_grid(plot.row.prop, legend.prop, rel_widths = c(2, .75))

# save pdf
# pdf(file = "Prop_cases_by_assumption.pdf", width = 7, height = 3.5, useDingbats=F)
# plot.fullprop
# dev.off()

# save eps for submission
save_plot("Fig2.eps", plot.fullprop,
          base_height = 3.5, base_width = 7)

# CHECK ADMISSION PREVALENCE  =======================================================

# find disease equilibrium
modpars <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
             gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
             rSH, rIL, rIH, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = 0, 
             gamma.intL, gamma.intH, decol = 0)
mod.init <- c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
              SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0, 
              SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
              SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0, Total.t0)
names(mod.init) <- c("SL1c", "SH1c", "IL1c", "IH1c", "DL1c", "DH1c",
                     "SL2c", "SH2c", "IL2c", "IH2c", "DL2c", "DH2c",
                     "SL1n", "SH1n", "IL1n", "IH1n", "DL1n", "DH1n",
                     "SL2n", "SH2n", "IL2n", "IH2n", "DL2n", "DH2n",
                     "total")
no.int.equil <- runsteady(y = mod.init, func = modsim.dyn, parms = modpars, times = c(0, 10000000))
equil.init <- unname(no.int.equil$y) # save like this so it's easier to use

eqprev.hosp.noint <- sum(no.int.equil$y[c(3, 4, 5, 6, 15, 16, 17, 18)])/N1.t0
eqprev.comm.noint <- sum(no.int.equil$y[c(9, 10, 11, 12, 21, 22, 23, 24)])/N2.t0

# simulate to get admission prevalence
modpars <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
             gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
             rSH, rIL, rIH, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = 0, 
             gamma.intL, gamma.intH, decol = 0)
mod.init <- c(SL1c.t0 = equil.init[1], SH1c.t0 = equil.init[2], IL1c.t0 = equil.init[3], 
              IH1c.t0 = equil.init[4], DL1c.t0 = equil.init[5], DH1c.t0 = equil.init[6], 
              SL2c.t0 = equil.init[7], SH2c.t0 = equil.init[8], IL2c.t0 = equil.init[9], 
              IH2c.t0 = equil.init[10], DL2c.t0 = equil.init[11], DH2c.t0 = equil.init[12], 
              SL1n.t0 = equil.init[13], SH1n.t0 = equil.init[14], IL1n.t0 = equil.init[15], 
              IH1n.t0 = equil.init[16], DL1n.t0 = equil.init[17], DH1n.t0 = equil.init[18], 
              SL2n.t0 = equil.init[19], SH2n.t0 = equil.init[20], IL2n.t0 = equil.init[21], 
              IH2n.t0 = equil.init[22], DL2n.t0 = equil.init[23], DH2n.t0 = equil.init[24],
              Total.t0 = equil.init[25])
mod.t <- seq(0, 1000, by=1) 
out <- as.data.frame(lsoda(mod.init, mod.t, modsim.dyn, modpars))
names(out) <- c("time", "SL1c", "SH1c", "IL1c", "IH1c", "DL1c", "DH1c",
                "SL2c", "SH2c", "IL2c", "IH2c", "DL2c", "DH2c",
                "SL1n", "SH1n", "IL1n", "IH1n", "DL1n", "DH1n",
                "SL2n", "SH2n", "IL2n", "IH2n", "DL2n", "DH2n", "Total",
                "Inc_IL1c", "Inc_IH1c", "Inc_IL2c", "Inc_IH2c", 
                "Inc_IL1n", "Inc_IH1n", "Inc_IL2n", "Inc_IH2n",
                "Inc_DL1c", "Inc_DH1c", "Inc_DL2c", "Inc_DH2c", 
                "Inc_DL1n", "Inc_DH1n", "Inc_DL2n", "Inc_DH2n",
                "adm_SLc", "adm_SLn", "adm_SHc", "adm_SHn", "adm_ILc", "adm_ILn",
                "adm_IHc", "adm_IHn", "adm_DLc", "adm_DLn", "adm_DHc", "adm_DHn",
                "IncS_IL1c", "IncS_IH1c", "IncS_IL1n", "IncS_IH1n")

# tail(out)

# number transmissions in hospital per 1,000 patient days
trans_hosp <- (sum(out[1000, c("Inc_IL1c", "Inc_IH1c", "Inc_IL1n", "Inc_IH1n")]))/N1.t0
trans1000days <- trans_hosp*1000

# number transmissions in hospital per 1,000 patient days (scaled by susceptible admissions)
trans_hosp_S <- (sum(out[1000, c("IncS_IL1c", "IncS_IH1c", "IncS_IL1n", "IncS_IH1n")]))
trans1000days_S <- trans_hosp_S*1000

# number of admissions by group
adm_ILc <- out[1000, "adm_ILc"]
adm_ILn <- out[1000, "adm_ILn"]
adm_IHc <- out[1000, "adm_IHc"]
adm_IHn <- out[1000, "adm_IHn"]
adm_DLc <- out[1000, "adm_DLc"]
adm_DLn <- out[1000, "adm_DLn"]
adm_DHc <- out[1000, "adm_DHc"]
adm_DHn <- out[1000, "adm_DHn"]

adm_Hc <- sum(out[1000, c("adm_SHc", "adm_IHc", "adm_DHc")])
adm_Hn <- sum(out[1000, c("adm_SHn", "adm_IHn", "adm_DHn")])
adm_Lc <- sum(out[1000, c("adm_SLc", "adm_ILc", "adm_DLc")])
adm_Ln <- sum(out[1000, c("adm_SLn", "adm_ILn", "adm_DLn")])

# proportion of admissions that are high risk (HCUP 2012 says this should be 0.349)
prop_adm_H <- (adm_Hc + adm_Hn)/(adm_Hc + adm_Hn + adm_Lc + adm_Ln)
prop_adm_H

# prevalence among newly admitted patients
adm_prev <- (adm_ILc + adm_ILn + adm_IHc + adm_IHn + adm_DLc + adm_DLn + adm_DHc + adm_DHn)/
  (adm_Hc + adm_Hn + adm_Lc + adm_Ln)

adm_H_prev <- (adm_IHc + adm_IHn + adm_DHc + adm_DHn)/(adm_Hc + adm_Hn)
adm_L_prev <- (adm_ILc + adm_ILn + adm_DLc + adm_DLn)/(adm_Lc + adm_Ln)


hosp.prev <- (out$IL1c + out$IH1c + out$DL1c + out$DH1c + out$IL1n + out$IH1n + out$DL1n + out$DH1n)/N1.t0
comm.prev <- (out$IL2c + out$IH2c + out$DL2c + out$DH2c + out$IL2n + out$IH2n + out$DL2n + out$DH2n)/N2.t0

out$hosp.prev <- hosp.prev
out$comm.prev <- comm.prev

max.yr <- round(max(out[,1])/365)
fifyearbreaks <- seq(0, max.yr, by = 50)
fifyearlabels <- seq(0, max.yr, by = 50)
text.spot <- round(max(out[,1])*0.75)

plot.transient.prev <- ggplot() +
  geom_line(data = out, aes(x = time, y = hosp.prev*100), color = "#7570b3", size = 1) +
  geom_line(data = out, aes(x = time, y = comm.prev*100), color = "#1b9e77", size = 1) +
  geom_hline(yintercept = 3.4, linetype = "dashed", color = "#7570b3", size = 1, alpha = 0.4) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "#1b9e77", size = 1, alpha = 0.4) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "#000000", size = 0.75, alpha = 0.4) +  
  ylab("prevalence (%)") +
  xlab("Time (years)") +
  ylim(0, 10) +
  scale_x_continuous(name = "Year", breaks = fifyearbreaks*365, labels = fifyearlabels) +
  annotate("text", x = text.spot, y = 9, label = paste("Admission prev = ", round(adm_prev*100, digits = 2),"%"))
plot.transient.prev


# GET WEIGHTED CLEARANCE RATES FOR HOMOGENEOUS ASSUMPTION  ==========================

# fit to parameters with no risk groups, but with long/short-term carriers
p <- 0.15
q <- 0.2
comm <- 550
old.beta1 <- 0.14102773
delta1 <- 0.34
beta2 <- 0.01096798
delta2 <- 0.167 
sigma <- 0.51
gamma1n <- 1/30
gamma1c <- 1/365
gamma2n <- 1/30
gamma2c <- 1/365
gammaD <- 1/5
alphaIL1 <- 0.009
alphaIH1 <- 0.009
alphaIL2 <- 0.009/5
alphaIH2 <- 0.009/5
rSL <- 1/4.5 # 34.8 may need to recalculate this based on admission rates of high and low risk patients, instead of p
rSH <- 1/4.5
rIL <- 1/4.5
rIH <- 1/4.5
rD <- 0.5 # multiplier
rhoH <- 1
rhoI <- 1
rhoD <- 2 # CHANGED
int <- 0 # not in sensitivity analysis
time.int <- 30000 # not in sensitivity analysis 
beta.int <- 0.30
gamma.intL <- 1/28
gamma.intH <- 1/28
decol <- 0 # no decoloniztion

# need to run the model with no infecteds to get proportions at equilibrium
N1.t0 <- 300
N2.t0 <- N1.t0*comm 
N1c.t0 <- q*N1.t0
N1n.t0 <- (1-q)*N1.t0
N2c.t0 <- q*N2.t0
N2n.t0 <- (1-q)*N2.t0

IL1c.t0 <- 0
IH1c.t0 <- 0
DL1c.t0 <- 0
DH1c.t0 <- 0
IL2c.t0 <- 0
IH2c.t0 <- 0
DL2c.t0 <- 0
DH2c.t0 <- 0
SL1c.t0 <- N1c.t0*(1-p) - IL1c.t0 - DL1c.t0
SH1c.t0 <- N1c.t0*(p) - IH1c.t0 - DH1c.t0
SL2c.t0 <- N2c.t0*(1-p) - IL2c.t0 - DL2c.t0
SH2c.t0 <- N2c.t0*(p) - IH2c.t0 - DH2c.t0

IL1n.t0 <- 0
IH1n.t0 <- 0
DL1n.t0 <- 0
DH1n.t0 <- 0
IL2n.t0 <- 0
IH2n.t0 <- 0
DL2n.t0 <- 0
DH2n.t0 <- 0
SL1n.t0 <- N1n.t0*(1-p) - IL1n.t0 - DL1n.t0
SH1n.t0 <- N1n.t0*(p) - IH1n.t0 - DH1n.t0
SL2n.t0 <- N2n.t0*(1-p) - IL2n.t0 - DL2n.t0
SH2n.t0 <- N2n.t0*(p) - IH2n.t0 - DH2n.t0

Total.t0 <- sum(c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
                  SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0,
                  SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
                  SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0))

# run to find equilibrium with no infecteds
modpars <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
             gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
             rSH, rIL, rIH, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = 0, 
             gamma.intL, gamma.intH, decol = 0)
mod.init <- c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
              SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0, 
              SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
              SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0, Total.t0)
names(mod.init) <- c("SL1c", "SH1c", "IL1c", "IH1c", "DL1c", "DH1c",
                     "SL2c", "SH2c", "IL2c", "IH2c", "DL2c", "DH2c",
                     "SL1n", "SH1n", "IL1n", "IH1n", "DL1n", "DH1n",
                     "SL2n", "SH2n", "IL2n", "IH2n", "DL2n", "DH2n",
                     "total")
uninf.equil <- runsteady(y = mod.init, func = modsim.dyn, parms = modpars, times = c(0, 10000000))

# size of hospital, community, etc
run_N1 <- sum(uninf.equil$y[c(1:6, 13:18)])
run_N2 <- sum(uninf.equil$y[c(7:12, 19:24)])

# size by risk group and location
run_NL1 <- sum(uninf.equil$y[c("SL1c", "IL1c", "DL1c", "SL1n", "IL1n", "DL1n")])
run_NH1 <- sum(uninf.equil$y[c("SH1c", "IH1c", "DH1c", "SH1n", "IH1n", "DH1n")])
run_NL2 <- sum(uninf.equil$y[c("SL2c", "IL2c", "DL2c", "SL2n", "IL2n", "DL2n")])
run_NH2 <- sum(uninf.equil$y[c("SH2c", "IH2c", "DH2c", "SH2n", "IH2n", "DH2n")])

# size by carrier type and location
run_N1c <- sum(uninf.equil$y[c("SL1c", "IL1c", "DL1c", "SH1c", "IH1c", "DH1c")])
run_N2c <- sum(uninf.equil$y[c("SL2c", "IL2c", "DL2c", "SH2c", "IH2c", "DH2c")])

# proportion high risk (p) and proportion recently discharged (q)
p1 <- run_NH1/run_N1
p2 <- run_NH2/run_N2
q1 <- run_N1c/run_N1
q2 <- run_N2c/run_N2


# fit to normal prev ================================================================

NL2.eq <- sum(uninf.equil$y[c(7, 9, 11, 19, 21, 23)])
NH2.eq <- sum(uninf.equil$y[c(8, 10, 12, 20, 22, 24)])
NL1.eq <- sum(uninf.equil$y[c(1, 3, 5, 13, 15, 17)])


N1.t0 <- 300
N2.t0 <- N1.t0*comm 
N1c.t0 <- q1*N1.t0
N1n.t0 <- (1-q1)*N1.t0
N2c.t0 <- q2*N2.t0
N2n.t0 <- (1-q2)*N2.t0

IL1c.t0 <- 1
IH1c.t0 <- 1
DL1c.t0 <- 0
DH1c.t0 <- 0
IL2c.t0 <- 1
IH2c.t0 <- 1
DL2c.t0 <- 0
DH2c.t0 <- 0
SL1c.t0 <- N1c.t0*(1-p1) - IL1c.t0 - DL1c.t0
SH1c.t0 <- N1c.t0*(p1) - IH1c.t0 - DH1c.t0
SL2c.t0 <- N2c.t0*(1-p2) - IL2c.t0 - DL2c.t0
SH2c.t0 <- N2c.t0*(p2) - IH2c.t0 - DH2c.t0

IL1n.t0 <- 1
IH1n.t0 <- 1
DL1n.t0 <- 0
DH1n.t0 <- 0
IL2n.t0 <- 1
IH2n.t0 <- 1
DL2n.t0 <- 0
DH2n.t0 <- 0
SL1n.t0 <- N1n.t0*(1-p1) - IL1n.t0 - DL1n.t0
SH1n.t0 <- N1n.t0*(p1) - IH1n.t0 - DH1n.t0
SL2n.t0 <- N2n.t0*(1-p2) - IL2n.t0 - DL2n.t0
SH2n.t0 <- N2n.t0*(p2) - IH2n.t0 - DH2n.t0

Total.t0 <- sum(c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
                  SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0,
                  SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
                  SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0))


# fit the transmission parameters
params0 <- c(log(0.1427436), log(0.01097243))  # initial guess, transformed so range is (0, +inf)
fit2 <- optim(params0, sse.rs.prev) # fit
beta.fit <- exp(fit2$par) # back transform

old.beta1 <- beta.fit[1]
beta2 <- beta.fit[2]

# equilibrium with no intervention
modpars <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
             gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
             rSH, rIL, rIH, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = 0, 
             gamma.intL, gamma.intH, decol = 0)
mod.init <- c(SL1c.t0, SH1c.t0, IL1c.t0, IH1c.t0, DL1c.t0, DH1c.t0, 
              SL2c.t0, SH2c.t0, IL2c.t0, IH2c.t0, DL2c.t0, DH2c.t0, 
              SL1n.t0, SH1n.t0, IL1n.t0, IH1n.t0, DL1n.t0, DH1n.t0, 
              SL2n.t0, SH2n.t0, IL2n.t0, IH2n.t0, DL2n.t0, DH2n.t0, Total.t0)
names(mod.init) <- c("SL1c", "SH1c", "IL1c", "IH1c", "DL1c", "DH1c",
                     "SL2c", "SH2c", "IL2c", "IH2c", "DL2c", "DH2c",
                     "SL1n", "SH1n", "IL1n", "IH1n", "DL1n", "DH1n",
                     "SL2n", "SH2n", "IL2n", "IH2n", "DL2n", "DH2n",
                     "total")
no.int.equil <- runsteady(y = mod.init, func = modsim.dyn, parms = modpars, times = c(0, 10000000))

eqprev.hosp.noint <- sum(no.int.equil$y[c(3, 4, 5, 6, 15, 16, 17, 18)])/N1.t0
eqprev.comm.noint <- sum(no.int.equil$y[c(9, 10, 11, 12, 21, 22, 23, 24)])/N2.t0
eqprev.hosp.noint; eqprev.comm.noint

# find proportions of colonized & diseased indivs that are in each carrier length group 
end.noint <- as.data.frame(t(no.int.equil$y))

eq_carriers <- no.int.equil$y[c(3:6, 9:12, 15:18, 21:24)]
eq_short_carriers <- eq_carriers[9:16]
eq_long_carriers <- eq_carriers[1:8]

sum_eq_carriars <- sum(eq_carriers) 
sum_eq_short_carriers <- sum(eq_short_carriers)
sum_eq_long_carriers <- sum(eq_long_carriers)

# fraction of carriers in short duration group
f_cs <- sum_eq_short_carriers/sum_eq_carriars
f_cl <- sum_eq_long_carriers/sum_eq_carriars

avg_carrier_carriage <- f_cs*(1/gamma1n) + f_cl*(1/gamma1c)
