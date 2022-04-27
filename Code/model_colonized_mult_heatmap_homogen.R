# comparing effect of each multiplier (LOS vs admission rate)
# *re-fitting beta1, beta2 for each new value

source("model_functions.R")
require(deSolve)
require(rootSolve)
require(tidyverse)
require(cowplot)
require(directlabels)

# set plot font size for everything
theme_set(theme_cowplot(font_size = 12))


# loop to fit beta1, beta2

# HOMOGENEOUS VERSION OF PARAMS
p <- 0.15
q <- 0.2
comm <- 550
old.beta1 <- 0.14102773
delta1 <- 0
beta2 <- 0.01096798
delta2 <- 0 
sigma <- 1
gamma1n <- 1/278
gamma1c <- 1/278
gamma2n <- 1/278
gamma2c <- 1/278
gammaD <- 1/5
alphaIL1 <- 0.009
alphaIH1 <- 0.009
alphaIL2 <- 0.009/5
alphaIH2 <- 0.009/5
rSL <- 1/4.5 
rH <- 1
rI <- 1
rD <- 0.5 
rhoH <- 1
rhoI <- 1
rhoD <- 2 
int <- 0 # not in sensitivity analysis
time.int <- 30000 # not in sensitivity analysis 
beta.int <- 0.30
gamma.intL <- 1/28
gamma.intH <- 1/28
decol <- 0 # not in sensitivity analysis 
# 
# homogen <- c(sigma = 1, delta1 = 0, delta2 = 0, rSL = 1/4.5, rH = 1, 
#              rI = 1, rhoH = 1, rhoI = 1, gamma1n = 1/278, 
#              gamma1c = 1/278, gamma2n = 1/278, gamma2c = 1/278)


# old fits of beta1, beta2
old_fit <- read.csv("fit_data3.csv")

#old_fit <- old_fit[,-1]

rhoI.vec <- seq(1, 3, length.out = 20)
rI.vec <- seq(1, 3, length.out = 20)

# rhoI.vec <- seq(1, 3, length.out = 20)
# rI.vec <- seq(1, 3, length.out = 20)

fit_data <- as.data.frame(matrix(rep(NA, length(rhoI.vec)*length(rI.vec)*6), nrow = length(rhoI.vec)*length(rI.vec)))
names(fit_data) <- c("rhoI", "rI", "beta1", "beta2", "hosp.noint", "comm.noint")

h <- 1

for (i in 1:length(rhoI.vec)){
  rhoI <- rhoI.vec[i]
  for (j in 1:length(rI.vec)){
    rI <- rI.vec[j]
    
    
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
    
    
    # if fitting to start:
    # fit the transmission parameters
    # params0 <- c(log(0.1), log(0.001))  # initial guess, transformed so range is (0, +inf)
    # fit2 <- optim(params0, sse.rs.prev) # fit
    # beta.fit <- exp(fit2$par) # back transform
    # 
    # old.beta1 <- beta.fit[1]
    # beta2 <- beta.fit[2]

    # if fitting based on old fits:
    # fit the transmission parameters
    params0 <- c(log(old_fit[h,"beta1"]), log(old_fit[h,"beta2"]))
    fit2 <- optim(params0, sse.rs.prev) # fit
    beta.fit <- exp(fit2$par) # back transform

    old.beta1 <- beta.fit[1]
    beta2 <- beta.fit[2]
    # 
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
    
    fit_data[h, 1] <- rhoI.vec[i]
    fit_data[h, 2] <- rI.vec[j]
    fit_data[h, 3] <- old.beta1
    fit_data[h, 4] <- beta2
    fit_data[h, 5] <- eqprev.hosp.noint
    fit_data[h, 6] <- eqprev.comm.noint
    
    print(fit_data[h, ])
    h <- h + 1
    
  }  
}


fit_data <- na.omit(fit_data)

write.csv(fit_data, file = "fit_data3.csv")

# fit_data <- read.csv("fit_data2.csv")

# if fitting is not w/in certain equil values, remove

# hosp.min <- 0.034 - 0.034*0.1 # 10% lower
# hosp.max <- 0.034 + 0.034*0.1 
# comm.min <- 0.015 - 0.015*0.1 # 10% lower
# comm.max <- 0.015 + 0.015*0.1 

hosp.min <- 0.0335
hosp.max <- 0.0345
comm.min <- 0.0145
comm.max <- 0.0155


str(refit_data)


fit_data <- fit_data %>% select(-X)

ttt <- fit_data %>%
  filter(hosp.noint <= hosp.max&hosp.noint >= hosp.min) %>%
  filter(comm.noint <= comm.max&comm.noint >= comm.min)
rrr <- refit_data %>%
  filter(hosp.noint <= hosp.max&hosp.noint >= hosp.min) %>%
  filter(comm.noint <= comm.max&comm.noint >= comm.min)

ttt <- rbind(ttt, rrr)



beta.data <- as.data.frame(matrix(rep(NA, 8*nrow(ttt)), nrow = nrow(ttt)))
names(beta.data) <- c("rhoI", "rI", "hospital", "community", "cases.averted", 
                      "baseline.cases", "intervention.cases", "prop.cases.averted")

j <- 1

# add transmission reduction
new.beta.int <- 0.3

for (i in 1:nrow(ttt)) {
  rhoI <- ttt[i, 1]
  rI <- ttt[i, 2]
  old.beta1 <- ttt[i, 3]
  beta2 <- ttt[i, 4]
  
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
  # new.beta1 <- old.beta1*(1 - new.beta.int)
  # R0etc <- calc_R0(comm, q1, q2, p1, p2, rSL, rH, rI, rD, rhoH, rhoI, rhoD, old.beta1 = new.beta1, int, delta1, beta2, sigma, delta2, 
  #                  gamma1c, alphaIL1, gamma1n, alphaIH1, gamma2c, alphaIL2, gamma2n, alphaIH2, gammaD)
  # 
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
  
  beta.data[j, 1] <- rhoI
  beta.data[j, 2] <- rI
  beta.data[j, 3] <- eqprev.hosp.int
  beta.data[j, 4] <- eqprev.comm.int
  beta.data[j, 5] <- cases.averted
  beta.data[j, 6] <- baseline.cases
  beta.data[j, 7] <- intervention.cases
  beta.data[j, 8] <- prop.cases.averted
  
  j <- j + 1
  print(j)
}


str(beta.data)

p.1.beta <- ggplot(beta.data,
                   aes(x = rhoI, y = rI, z = cases.averted, color = cases.averted)) +
  geom_raster(aes(fill = cases.averted), interpolate = FALSE) +
  scale_fill_viridis_c() +
  stat_contour(bins = 10, color = "white") 
p.1.beta

p.2.beta <- ggplot(beta.data,
                   aes(x = rhoI, y = rI, z = baseline.cases, color = baseline.cases)) +
  geom_raster(aes(fill = baseline.cases), interpolate = FALSE) +
  scale_fill_viridis_c() +
  stat_contour(bins = 10, color = "white") 
p.2.beta

p.3.beta <- ggplot(beta.data,
                   aes(x = rhoI, y = rI, z = prop.cases.averted*100, color = prop.cases.averted*100)) +
  geom_raster(aes(fill = prop.cases.averted*100), interpolate = FALSE) +
  scale_fill_viridis_c() +
  stat_contour(bins = 10, color = "white") 
p.3.beta


# adjusting data sets to compare intervention to multiplier=1 scenario

str(beta.data)

mult1.noint.cases <- beta.data[1,"baseline.cases"]
mult1.int.cases <- beta.data[1,"intervention.cases"]

# add columns for cases when multipliers = 1 (with and w/o intervention)
beta.data$mult1.noint.cases <- rep(mult1.noint.cases, length.out = nrow(beta.data))
beta.data$mult1.int.cases <- rep(mult1.int.cases, length.out = nrow(beta.data))


test <- beta.data %>%
  mutate(rel.value = 100*(intervention.cases/baseline.cases - mult1.int.cases/mult1.noint.cases))

p.4.beta <- ggplot(test,
                   aes(x = rhoI, y = rI, z = rel.value, color = rel.value)) +
  geom_raster(aes(fill = rel.value), interpolate = FALSE) +
  scale_fill_viridis_c() +
  stat_contour(bins = 10, color = "white") +
  labs(fill = "Percent change\ncases averted")
#  theme(legend.title = element_text("Percent change\ncases averted"))
p.4.beta



# decolonization loop
decol.data <- as.data.frame(matrix(rep(NA, 8*nrow(ttt)), nrow = nrow(ttt)))
names(decol.data) <- c("rhoI", "rI", "hospital", "community", "cases.averted", 
                       "baseline.cases", "intervention.cases", "prop.cases.averted")

j <- 1

# add decolonization reduction
gamma.intH <- 0.035
gamma.intL <- 0.035

for (i in 1:nrow(ttt)) {
  rhoI <- ttt[i, 1]
  rI <- ttt[i, 2]
  old.beta1 <- ttt[i, 3]
  beta2 <- ttt[i, 4]
  
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
  # new.beta1 <- old.beta1*(1 - new.beta.int)
  # R0etc <- calc_R0(comm, q1, q2, p1, p2, rSL, rH, rI, rD, rhoH, rhoI, rhoD, old.beta1 = new.beta1, int, delta1, beta2, sigma, delta2, 
  #                  gamma1c, alphaIL1, gamma1n, alphaIH1, gamma2c, alphaIL2, gamma2n, alphaIH2, gammaD)
  # 
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
  
  decol.data[j, 1] <- rhoI
  decol.data[j, 2] <- rI
  decol.data[j, 3] <- eqprev.hosp.int
  decol.data[j, 4] <- eqprev.comm.int
  decol.data[j, 5] <- cases.averted
  decol.data[j, 6] <- baseline.cases
  decol.data[j, 7] <- intervention.cases
  decol.data[j, 8] <- prop.cases.averted
  
  j <- j + 1
  print(j)
}

p.1.decol <- ggplot(decol.data,
                    aes(x = rhoI, y = rI, z = cases.averted, color = cases.averted)) +
  geom_raster(aes(fill = cases.averted), interpolate = FALSE) +
  scale_fill_viridis_c() +
  stat_contour(bins = 10, color = "white") 
p.1.decol

p.2.decol <- ggplot(decol.data,
                    aes(x = rhoI, y = rI, z = baseline.cases, color = baseline.cases)) +
  geom_raster(aes(fill = baseline.cases), interpolate = FALSE) +
  scale_fill_viridis_c() +
  stat_contour(bins = 10, color = "white") 
p.2.decol

p.3.decol <- ggplot(decol.data,
                    aes(x = rhoI, y = rI, z = prop.cases.averted*100, color = prop.cases.averted*100)) +
  geom_raster(aes(fill = prop.cases.averted*100), interpolate = FALSE) +
  scale_fill_viridis_c() +
  stat_contour(bins = 10, color = "white") +
  labs(fill = "Percent change\ncases averted")
#  theme(legend.title = element_text("Percent change\ncases averted"))
p.3.decol



# adjusting data sets to compare intervention to multiplier=1 scenario

str(decol.data)

mult1.noint.cases <- decol.data[1,"baseline.cases"]
mult1.int.cases <- decol.data[1,"intervention.cases"]

# add columns for cases when multipliers = 1 (with and w/o intervention)
decol.data$mult1.noint.cases <- rep(mult1.noint.cases, length.out = nrow(decol.data))
decol.data$mult1.int.cases <- rep(mult1.int.cases, length.out = nrow(decol.data))


test <- decol.data %>%
  mutate(rel.value = 100*(intervention.cases/baseline.cases - mult1.int.cases/mult1.noint.cases))

p.4.decol <- ggplot(test,
                    aes(x = rhoI, y = rI, z = rel.value, color = rel.value)) +
  geom_raster(aes(fill = rel.value), interpolate = FALSE) +
  scale_fill_viridis_c() +
  stat_contour(bins = 10, color = "white") +
  labs(fill = "Percent change\ncases averted")
#  theme(legend.title = element_text("Percent change\ncases averted"))
p.4.decol



beta.data$treatment <- rep("transmission reduction", length.out = nrow(beta.data))
decol.data$treatment <- rep("decolonization", length.out = nrow(decol.data))

str(beta.data)
str(decol.data)

intervention.data <- bind_rows(beta.data, decol.data)

intervention.data <- intervention.data %>%
  mutate(rel.value = 100*(intervention.cases/baseline.cases - mult1.int.cases/mult1.noint.cases))

intervention.data$treatment <- as.factor(intervention.data$treatment)

write.csv(intervention.data, file = "data_multiplier_casesaverted.csv", row.names = FALSE)

intervention.data <- read.csv("data_multiplier_casesaverted.csv")


# new risk ratio metric =====================================================
str(intervention.data)


intervention.RR <- intervention.data %>%
  mutate(RR.value = 100*(mult1.int.cases/mult1.noint.cases - intervention.cases/baseline.cases)/
           (1 - mult1.int.cases/mult1.noint.cases))

summary(intervention.RR)


# making this so we can fill in NA's for combinations where fitting didn't work
# combos <-expand.grid(rhoI.vec, rI.vec)
# colnames(combos) <- c("rhoI", "rI")
# all.data <- full_join(intervention.RR, combos)


p.2.basic <- ggplot(intervention.RR,
                    aes(x = rhoI, y = rI,  z = RR.value)) +
  geom_tile(aes(fill = RR.value)) +
  scale_fill_viridis_c(direction = -1) +
  labs(fill = "% change in\ncases averted") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab(bquote("Colonized LOS multiplier," ~~ r[C])) +
  xlab(bquote("Colonized admission rate multiplier," ~~ rho[C])) +
  panel_border(size = 1, colour = "black") +
  theme(aspect.ratio = 1) +
  facet_wrap(~treatment) 
p.2.basic

p.2.both <- ggplot(intervention.RR,
                   aes(x = rhoI, y = rI,  z = RR.value)) +
  geom_tile(aes(fill = RR.value)) +
  scale_fill_viridis_c(direction = -1) +
  scale_colour_gradientn(colours = c("grey","grey"),
                         values = c(0,100)) +
  # stat_contour(bins = 10, color = "white") +
  stat_contour(aes(colour = ..level..)) +
  labs(fill = "% change in\ncases averted") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab(bquote("Colonized LOS multiplier," ~~ r[C])) +
  xlab(bquote("Colonized admission rate multiplier," ~~ rho[C])) +
  panel_border(size = 1, colour = "black") +
  theme(aspect.ratio = 1, legend.position = "none") +
  facet_wrap(~treatment) 
#p.2.both




p.2.leg <- get_legend(p.2.basic)


p.2.withlabs <- direct.label(p.2.both, list("bottom.pieces",
                                            hjust = 0.5, 
                                            vjust = 0.8,
                                            cex = 0.75, color = "white"))


full<-plot_grid(p.2.withlabs, p.2.leg, ncol = 2,
                rel_widths = c(2, .3),
                rel_heights = c(2, .3))
full

# pdf(file = "RR_homogen_Multipliers_cases_intervention.pdf", width = 8, height = 4, useDingbats=F)
# full
# dev.off()

# save eps for submission
save_plot("S2_Fig.eps", fullheatmap,
          base_height = 4, base_width = 8)



# "raster" version
p.1.both <- ggplot(intervention.data,
                   aes(x = rhoI, y = rI, z = rel.value, color = rel.value)) +
  geom_raster(aes(fill = rel.value), interpolate = FALSE) +
  scale_fill_viridis_c() +
  stat_contour(bins = 10, color = "white") +
  labs(fill = "Percent change\ncases averted") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab(bquote("Colonized LOS multiplier," ~~ r[C])) +
  xlab(bquote("Colonized admission rate multiplier," ~~ rho[C])) +
  panel_border(size = 1, colour = "black") +
  theme(aspect.ratio = 1) +
  facet_wrap(~treatment) 
p.1.both

# "tile" version
p.2.both <- ggplot(intervention.data,
                   aes(x = rhoI, y = rI,  z = rel.value)) +
  geom_tile(aes(fill = rel.value), color = NA) +
  scale_fill_viridis_c() +
  stat_contour(bins = 10, color = "white") +
  labs(fill = "Percent change\ncases averted") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab(bquote("Colonized LOS multiplier," ~~ r[C])) +
  xlab(bquote("Colonized admission rate multiplier," ~~ rho[C])) +
  panel_border(size = 1, colour = "black") +
  theme(aspect.ratio = 1) +
  facet_wrap(~treatment) 
p.2.both


# save plots as pdfs
# pdf(file = "Multipliers_cases_intervention.pdf", width = 8, height = 6, useDingbats=F)
# p.2.both
# dev.off()

str(intervention.data)

levels(intervention.data$treatment)

# individual plots
p.2.trans <- ggplot(intervention.data %>% filter(treatment == "transmission reduction"),
                    aes(x = rhoI, y = rI,  z = rel.value)) +
  geom_tile(aes(fill = rel.value), color = NA) +
  scale_fill_viridis_c() +
  stat_contour(bins = 10, color = "white") +
  labs(fill = "Percent change\ncases averted") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab(bquote("Colonized LOS multiplier," ~~ r[C])) +
  xlab(bquote("Colonized admission rate multiplier," ~~ rho[C])) +
  panel_border(size = 1, colour = "black") +
  ggtitle("transmission reduction") +
  theme(aspect.ratio = 1) 
p.2.trans

p.2.decol <- ggplot(intervention.data %>% filter(treatment == "decolonization"),
                    aes(x = rhoI, y = rI,  z = rel.value)) +
  geom_tile(aes(fill = rel.value), color = NA) +
  scale_fill_viridis_c() +
  stat_contour(bins = 10, color = "white") +
  labs(fill = "Percent change\ncases averted") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab(bquote("Colonized LOS multiplier," ~~ r[C])) +
  xlab(bquote("Colonized admission rate multiplier," ~~ rho[C])) +
  panel_border(size = 1, colour = "black") +
  ggtitle("transmission reduction") +
  theme(aspect.ratio = 1) 
p.2.decol


# save plots as pdfs
# pdf(file = "Multipliers_cases_intervention_sep.pdf", width = 6, height = 6, useDingbats=F)
# p.2.trans
# p.2.decol
# dev.off()


# loop to try re-fitting betas where it failed before
# take beta values from similar rH/rI values and start fitting from there



to_refit <- fit_data %>%
  mutate(refits = case_when(hosp.noint <= hosp.max&hosp.noint >= hosp.min ~ "N",
                            comm.noint <= comm.max&comm.noint >= comm.min ~ "N"))
to_refit[is.na(to_refit)] <- "Y"

to_refit <- to_refit %>% select(-X)

# read in old data

refit_data <- read.csv("refit_data2.csv")
refit_data <- na.omit(refit_data)

str(ttt)

str(refit_data)

hist(refit_data$hosp.noint)
abline(v = hosp.min, col = "red")
abline(v = hosp.max, col = "red")

hist(refit_data$comm.noint)
abline(v = comm.min, col = "red")
abline(v = comm.max, col = "red")


# rH.vec <- seq(1, 3, length.out = 10)
# rI.vec <- seq(1, 3, length.out = 10)

refit_data <- as.data.frame(matrix(rep(NA, length(rhoI.vec)*length(rI.vec)*6), nrow = length(rhoI.vec)*length(rI.vec)))
names(refit_data) <- c("rhoI", "rI", "beta1", "beta2", "hosp.noint", "comm.noint")

h <- 1

for (i in 1:nrow(to_refit)){
  
  if (to_refit[i, 7] == "Y") {
    rhoI <- to_refit[i, 1]
    rI <- to_refit[i, 2]
    
    # take params from previous value
    params0 <- c(log(to_refit[i-1, 3]), log(to_refit[i-1, 4]))  
    
    
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
    # params0 <- c(log(0.1), log(0.01))  # initial guess, transformed so range is (0, +inf)
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
    
    refit_data[h, 1] <- rhoI
    refit_data[h, 2] <- rI
    refit_data[h, 3] <- old.beta1
    refit_data[h, 4] <- beta2
    refit_data[h, 5] <- eqprev.hosp.noint
    refit_data[h, 6] <- eqprev.comm.noint
    
    h <- h + 1
    print(h)
    
  }
}



