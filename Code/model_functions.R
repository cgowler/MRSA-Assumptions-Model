# ODE Model, R0 calculation, and Parameter fitting equations

# Variable naming conventions are different than appear in the manuscript:
# "S" = "susceptible" (model) --> "U" = "uncolonized" (manuscript)
# "I" = "infectious" (model) --> "C" = "colonized" (manuscript)
# "D" = "diseased" (model) --> "S" = "symptomatic" (manuscript)

# MODEL  ==========================================================================
modsim.dyn <- function(t, var, par) { 
  
  # Variables
  SL1c <- var[1] 
  SH1c <- var[2] 
  IL1c <- var[3] 
  IH1c <- var[4] 
  DL1c <- var[5] 
  DH1c <- var[6] 
  SL2c <- var[7] 
  SH2c <- var[8] 
  IL2c <- var[9] 
  IH2c <- var[10] 
  DL2c <- var[11] 
  DH2c <- var[12] 
  SL1n <- var[13] 
  SH1n <- var[14] 
  IL1n <- var[15] 
  IH1n <- var[16] 
  DL1n <- var[17] 
  DH1n <- var[18] 
  SL2n <- var[19] 
  SH2n <- var[20]
  IL2n <- var[21] 
  IH2n <- var[22]
  DL2n <- var[23] 
  DH2n <- var[24]
  Total <- var[25]
  
  # Parameters
  old.beta1 <- par[1] # transmission in hospital (prior to intervention)
  delta1<- par[2] # hospital assortativity (1 = assortative, 0 = random)
  beta2 <- par[3] # transmission in community
  delta2 <- par[4] # community assortativity (1 = assortative, 0 = random)
  sigma <- par[5] # reduction in contacts of high risk group compared to low risk group in community
  gamma1n <- par[6] # recovery of non-carrier (Hospital)
  gamma1c <- par[7] # recovery of carrier (Hospital)
  gamma2n <- par[8] # recovery of non-carrier (Community)
  gamma2c <- par[9] # recovery of carrier (Community)
  gammaD   <- par[10] # disease state recovery rate 
  alphaIL1 <- par[11] # disease progression rate for low-risk in hospital
  alphaIH1 <- par[12] # disease progression rate for high-risk in hospital
  alphaIL2 <- par[13] # disease progression rate for low-risk in community
  alphaIH2 <- par[14] # disease progression rate for high-risk in community
  rSL <- par[15] # discharge rate of low-risk susceptible
  rH <- par[16] # discharge rate multiplier for high-risk individuals
  rI <- par[17] # discharge rate multiplier for infectious individuals
  rD <- par[18] # discharge rate multiplier for disease individuals
  rhoH <- par[19] # multiplier for increase in hospitalization rate for high risk 
  rhoI <- par[20] # multiplier for increase in hospitalization rate for infectious state
  rhoD <- par[21] # multiplier for increase in hospitalization rate for disease state
  int <- par[22] # intervention efficacy prior to intervention (always set to 0)
  time.int <- par[23] # time of intervention
  beta.int <- par[24] # intervention efficacy for transmission reduction after intervention
  gamma.intL <- par[25] # new recovery rate for decolonization treatment 
  gamma.intH <- par[26] # new recovery rate for decolonization treatment 
  decol <- par[27] # 0 or 1 to turn off/on decolonization treatment
  
  # Intervention conditions
  if(t >= time.int) int <- beta.int 
  if(t >= time.int & decol == 1) {gamma1Lc <- gamma1c + gamma.intL} else {gamma1Lc <- gamma1c} 
  if(t >= time.int & decol == 1) {gamma1Hc <- gamma1c + gamma.intH} else {gamma1Hc <- gamma1c} 
  if(t >= time.int & decol == 1) {gamma1Ln <- gamma1n + gamma.intL} else {gamma1Ln <- gamma1n} 
  if(t >= time.int & decol == 1) {gamma1Hn <- gamma1n + gamma.intH} else {gamma1Hn <- gamma1n}
  beta1 <- old.beta1*(1 - int)
  
  # Total population size 
  N1  <- sum(c(SL1c, SH1c, IL1c, IH1c, DL1c, DH1c, SL1n, SH1n, IL1n, IH1n, DL1n, DH1n))
  N2  <- sum(c(SL2c, SH2c, IL2c, IH2c, DL2c, DH2c, SL2n, SH2n, IL2n, IH2n, DL2n, DH2n))
  NL1 <- sum(c(SL1c, IL1c, DL1c, SL1n, IL1n, DL1n))
  NH1 <- sum(c(SH1c, IH1c, DH1c, SH1n, IH1n, DH1n))
  NL2 <- sum(c(SL2c, IL2c, DL2c, SL2n, IL2n, DL2n))
  NH2 <- sum(c(SH2c, IH2c, DH2c, SH2n, IH2n, DH2n))
  
  # Transmission groups
  ID_L1 <- IL1c + IL1n + DL1c + DL1n
  ID_H1 <- IH1c + IH1n + DH1c + DH1n
  ID_L2 <- IL2c + IL2n + DL2c + DL2n
  ID_H2 <- IH2c + IH2n + DH2c + DH2n
  
  # Balancing total discharge and admission rates
  rSH <- rSL/rH # need to do it this way because multiplier is meant to multiply the LOS not the rate!
  rIL <- rSL/rI
  rIH <- rSL/(rI*rH)
  rDL <- rD*rIL # this multiplier is meant for the rate so it's okay
  rDH <- rD*rIH
  
  hosp.discharge.rate <- (rSL*SL1c + rSH*SH1c + rIL*IL1c + rIH*IH1c + rDL*DL1c + rDH*DH1c +
                            rSL*SL1n + rSH*SH1n + rIL*IL1n + rIH*IH1n + rDL*DL1n + rDH*DH1n) 
  
  hospitalization.rate <- (hosp.discharge.rate)/(SL2c + rhoH*SH2c + rhoI*IL2c  + rhoH*rhoI*IH2c + rhoD*rhoI*DL2c + rhoD*rhoH*rhoI*DH2c +
                                                   SL2n + rhoH*SH2n + rhoI*IL2n  + rhoH*rhoI*IH2n + rhoD*rhoI*DL2n + rhoD*rhoH*rhoI*DH2n)
  aSL <- hospitalization.rate 
  aSH <- rhoH*hospitalization.rate
  aIL <- rhoI*hospitalization.rate
  aIH <- rhoH*rhoI*hospitalization.rate
  aDL <- rhoD*rhoI*hospitalization.rate
  aDH <- rhoD*rhoH*rhoI*hospitalization.rate
  
  # ODEs 
  dSL1c <- - beta1*((1 - delta1)*NH1/N1)*SL1c*(ID_H1/NH1) - beta1*(delta1 + (1 - delta1)*NL1/N1)*SL1c*(ID_L1/NL1) + gamma1Lc*IL1c + aSL*SL2c - rSL*SL1c
  dIL1c <- + beta1*((1 - delta1)*NH1/N1)*SL1c*(ID_H1/NH1) + beta1*(delta1 + (1 - delta1)*NL1/N1)*SL1c*(ID_L1/NL1) - gamma1Lc*IL1c + aIL*IL2c - rIL*IL1c - alphaIL1*IL1c + gammaD*DL1c
  dDL1c <- + alphaIL1*IL1c - gammaD*DL1c + aDL*DL2c - rDL*DL1c
  
  dSL1n <- - beta1*((1 - delta1)*NH1/N1)*SL1n*(ID_H1/NH1) - beta1*(delta1 + (1 - delta1)*NL1/N1)*SL1n*(ID_L1/NL1) + gamma1Ln*IL1n + aSL*SL2n - rSL*SL1n
  dIL1n <- + beta1*((1 - delta1)*NH1/N1)*SL1n*(ID_H1/NH1) + beta1*(delta1 + (1 - delta1)*NL1/N1)*SL1n*(ID_L1/NL1) - gamma1Ln*IL1n + aIL*IL2n - rIL*IL1n - alphaIL1*IL1n + gammaD*DL1n
  dDL1n <- + alphaIL1*IL1n - gammaD*DL1n + aDL*DL2n - rDL*DL1n
  
  dSH1c <- - beta1*((1 - delta1)*NL1/N1)*SH1c*(ID_L1/NL1) - beta1*(delta1 + (1 - delta1)*NH1/N1)*SH1c*(ID_H1/NH1) + gamma1Hc*IH1c + aSH*SH2c - rSH*SH1c
  dIH1c <- + beta1*((1 - delta1)*NL1/N1)*SH1c*(ID_L1/NL1) + beta1*(delta1 + (1 - delta1)*NH1/N1)*SH1c*(ID_H1/NH1) - gamma1Hc*IH1c + aIH*IH2c - rIH*IH1c - alphaIH1*IH1c + gammaD*DH1c
  dDH1c <- + alphaIH1*IH1c - gammaD*DH1c + aDH*DH2c - rDH*DH1c
  
  dSH1n <- - beta1*((1 - delta1)*NL1/N1)*SH1n*(ID_L1/NL1) - beta1*(delta1 + (1 - delta1)*NH1/N1)*SH1n*(ID_H1/NH1) + gamma1Hn*IH1n + aSH*SH2n - rSH*SH1n
  dIH1n <- + beta1*((1 - delta1)*NL1/N1)*SH1n*(ID_L1/NL1) + beta1*(delta1 + (1 - delta1)*NH1/N1)*SH1n*(ID_H1/NH1) - gamma1Hn*IH1n + aIH*IH2n - rIH*IH1n - alphaIH1*IH1n + gammaD*DH1n
  dDH1n <- + alphaIH1*IH1n - gammaD*DH1n + aDH*DH2n - rDH*DH1n
  
  dSL2c <- - beta2*sigma*((1 - delta2)*NH2/N2)*SL2c*(ID_H2/NH2) - beta2*(delta2 + (1 - delta2)*NL2/N2)*SL2c*(ID_L2/NL2) + gamma2c*IL2c - aSL*SL2c + rSL*SL1c
  dIL2c <- + beta2*sigma*((1 - delta2)*NH2/N2)*SL2c*(ID_H2/NH2) + beta2*(delta2 + (1 - delta2)*NL2/N2)*SL2c*(ID_L2/NL2) - gamma2c*IL2c - aIL*IL2c + rIL*IL1c - alphaIL2*IL2c + gammaD*DL2c
  dDL2c <- + alphaIL2*IL2c - gammaD*DL2c - aDL*DL2c + rDL*DL1c
  
  dSL2n <- - beta2*sigma*((1 - delta2)*NH2/N2)*SL2n*(ID_H2/NH2) - beta2*(delta2 + (1 - delta2)*NL2/N2)*SL2n*(ID_L2/NL2) + gamma2n*IL2n - aSL*SL2n + rSL*SL1n
  dIL2n <- + beta2*sigma*((1 - delta2)*NH2/N2)*SL2n*(ID_H2/NH2) + beta2*(delta2 + (1 - delta2)*NL2/N2)*SL2n*(ID_L2/NL2) - gamma2n*IL2n - aIL*IL2n + rIL*IL1n - alphaIL2*IL2n + gammaD*DL2n
  dDL2n <- + alphaIL2*IL2n - gammaD*DL2n - aDL*DL2n + rDL*DL1n
  
  dSH2c <- - beta2*((1 - delta2)*NL2/N2)*SH2c*(ID_L2/NL2) - beta2*sigma*(delta2 + (1 - delta2)*NH2/N2)*SH2c*(ID_H2/NH2) + gamma2c*IH2c - aSH*SH2c + rSH*SH1c
  dIH2c <- + beta2*((1 - delta2)*NL2/N2)*SH2c*(ID_L2/NL2) + beta2*sigma*(delta2 + (1 - delta2)*NH2/N2)*SH2c*(ID_H2/NH2) - gamma2c*IH2c - aIH*IH2c + rIH*IH1c - alphaIH2*IH2c + gammaD*DH2c
  dDH2c <- + alphaIH2*IH2c - gammaD*DH2c - aDH*DH2c + rDH*DH1c
  
  dSH2n <- - beta2*((1 - delta2)*NL2/N2)*SH2n*(ID_L2/NL2) - beta2*sigma*(delta2 + (1 - delta2)*NH2/N2)*SH2n*(ID_H2/NH2) + gamma2n*IH2n - aSH*SH2n + rSH*SH1n
  dIH2n <- + beta2*((1 - delta2)*NL2/N2)*SH2n*(ID_L2/NL2) + beta2*sigma*(delta2 + (1 - delta2)*NH2/N2)*SH2n*(ID_H2/NH2) - gamma2n*IH2n - aIH*IH2n + rIH*IH1n - alphaIH2*IH2n + gammaD*DH2n
  dDH2n <- + alphaIH2*IH2n - gammaD*DH2n - aDH*DH2n + rDH*DH1n
  
  dTotal <- (dSL1c + dSH1c + dIL1c + dIH1c + dSL2c + dSH2c + dIL2c + dIH2c + dDL1c + dDH1c + dDL2c + dDH2c +
               dSL1n + dSH1n + dIL1n + dIH1n + dSL2n + dSH2n + dIL2n + dIH2n + dDL1n + dDH1n + dDL2n + dDH2n)
  
  # Admission numbers for different groups 
  adm_SLc <- aSL*SL2c
  adm_SLn <- aSL*SL2n
  adm_SHc <- aSH*SH2c
  adm_SHn <- aSH*SH2n
  adm_ILc <- aIL*IL2c
  adm_ILn <- aIL*IL2n
  adm_IHc <- aIH*IH2c
  adm_IHn <- aIH*IH2n
  adm_DLc <- aDL*DL2c
  adm_DLn <- aDL*DL2n
  adm_DHc <- aDH*DH2c
  adm_DHn <- aDH*DH2n
  adm_S <- adm_SLc + adm_SLn + adm_SHc + adm_SHn
  
  # Incidence of infection 
  Inc_IL1c <- + beta1*((1 - delta1)*NH1/N1)*SL1c*(ID_H1/NH1) + beta1*(delta1 + (1 - delta1)*NL1/N1)*SL1c*(ID_L1/NL1)
  Inc_IH1c <- + beta1*((1 - delta1)*NL1/N1)*SH1c*(ID_L1/NL1) + beta1*(delta1 + (1 - delta1)*NH1/N1)*SH1c*(ID_H1/NH1)
  Inc_IL2c <- + beta2*sigma*((1 - delta2)*NH2/N2)*SL2c*(ID_H2/NH2) + beta2*(delta2 + (1 - delta2)*NL2/N2)*SL2c*(ID_L2/NL2)
  Inc_IH2c <- + beta2*((1 - delta2)*NL2/N2)*SH2c*(ID_L2/NL2) + beta2*sigma*(delta2 + (1 - delta2)*NH2/N2)*SH2c*(ID_H2/NH2)
  Inc_IL1n <- + beta1*((1 - delta1)*NH1/N1)*SL1n*(ID_H1/NH1) + beta1*(delta1 + (1 - delta1)*NL1/N1)*SL1n*(ID_L1/NL1)
  Inc_IH1n <- + beta1*((1 - delta1)*NL1/N1)*SH1n*(ID_L1/NL1) + beta1*(delta1 + (1 - delta1)*NH1/N1)*SH1n*(ID_H1/NH1)
  Inc_IL2n <- + beta2*sigma*((1 - delta2)*NH2/N2)*SL2n*(ID_H2/NH2) + beta2*(delta2 + (1 - delta2)*NL2/N2)*SL2n*(ID_L2/NL2)
  Inc_IH2n <- + beta2*((1 - delta2)*NL2/N2)*SH2n*(ID_L2/NL2) + beta2*sigma*(delta2 + (1 - delta2)*NH2/N2)*SH2n*(ID_H2/NH2)
  
  # Incidence of infection (scale by susceptibles only)
  IncS_IL1c <- (beta1*((1 - delta1)*NH1/N1)*adm_SLc*(ID_H1/NH1) + beta1*(delta1 + (1 - delta1)*NL1/N1)*adm_SLc*(ID_L1/NL1))/adm_S
  IncS_IH1c <- (beta1*((1 - delta1)*NL1/N1)*adm_SHc*(ID_L1/NL1) + beta1*(delta1 + (1 - delta1)*NH1/N1)*adm_SHc*(ID_H1/NH1))/adm_S
  IncS_IL1n <- (beta1*((1 - delta1)*NH1/N1)*adm_SLn*(ID_H1/NH1) + beta1*(delta1 + (1 - delta1)*NL1/N1)*adm_SLn*(ID_L1/NL1))/adm_S
  IncS_IH1n <- (beta1*((1 - delta1)*NL1/N1)*adm_SHn*(ID_L1/NL1) + beta1*(delta1 + (1 - delta1)*NH1/N1)*adm_SHn*(ID_H1/NH1))/adm_S
  
  # Incidence of disease
  Inc_DL1c <- alphaIL1*IL1c 
  Inc_DH1c <- alphaIH1*IH1c 
  Inc_DL2c <- alphaIL2*IL2c  
  Inc_DH2c <- alphaIH2*IH2c
  Inc_DL1n <- alphaIL1*IL1n 
  Inc_DH1n <- alphaIH1*IH1n 
  Inc_DL2n <- alphaIL2*IL2n  
  Inc_DH2n <- alphaIH2*IH2n
  
  list(c(dSL1c, dSH1c, dIL1c, dIH1c, dDL1c, dDH1c,
         dSL2c, dSH2c, dIL2c, dIH2c, dDL2c, dDH2c,
         dSL1n, dSH1n, dIL1n, dIH1n, dDL1n, dDH1n,
         dSL2n, dSH2n, dIL2n, dIH2n, dDL2n, dDH2n, dTotal),
       Inc_IL1c, Inc_IH1c, Inc_IL2c, Inc_IH2c, 
       Inc_IL1n, Inc_IH1n, Inc_IL2n, Inc_IH2n,
       Inc_DL1c, Inc_DH1c, Inc_DL2c, Inc_DH2c, 
       Inc_DL1n, Inc_DH1n, Inc_DL2n, Inc_DH2n,
       adm_SLc, adm_SLn, adm_SHc, adm_SHn, adm_ILc, adm_ILn,
       adm_IHc, adm_IHn, adm_DLc, adm_DLn, adm_DHc, adm_DHn, 
       IncS_IL1c, IncS_IH1c, IncS_IL1n, IncS_IH1n,
       gamma1Lc, gamma1Hc, gamma1Ln, gamma1Hn)
}

# R0  =============================================================================
calc_R0 <- function(comm, q1, q2, p1, p2, rSL, rH, rI, rD, rhoH, rhoI, rhoD, 
                    old.beta1, int, delta1, beta2, sigma, delta2, gamma1c, alphaIL1, 
                    gamma1n, alphaIH1, gamma2c, alphaIL2, gamma2n, alphaIH2, gammaD) {
  
  pop1 <- 300
  pop2 <- pop1*comm 
  pop1c <- q1*pop1
  pop1n <- (1-q1)*pop1
  pop2c <- q2*pop2
  pop2n <- (1-q2)*pop2
  
  IL1c.init <- 0
  IH1c.init <- 0
  DL1c.init <- 0
  DH1c.init <- 0
  IL2c.init <- 0
  IH2c.init <- 0
  DL2c.init <- 0
  DH2c.init <- 0
  SL1c.init <- pop1c*(1-p1) - IL1c.init - DL1c.init
  SH1c.init <- pop1c*(p1) - IH1c.init - DH1c.init
  SL2c.init <- pop2c*(1-p2) - IL2c.init - DL2c.init
  SH2c.init <- pop2c*(p2) - IH2c.init - DH2c.init
  
  IL1n.init <- 0
  IH1n.init <- 0
  DL1n.init <- 0
  DH1n.init <- 0
  IL2n.init <- 0
  IH2n.init <- 0
  DL2n.init <- 0
  DH2n.init <- 0
  SL1n.init <- pop1n*(1-p1) - IL1n.init - DL1n.init
  SH1n.init <- pop1n*(p1) - IH1n.init - DH1n.init
  SL2n.init <- pop2n*(1-p2) - IL2n.init - DL2n.init
  SH2n.init <- pop2n*(p2) - IH2n.init - DH2n.init
  
  popL1 <- sum(c(SL1c.init, IL1c.init, DL1c.init, SL1n.init, IL1n.init, DL1n.init))  
  popH1 <- sum(c(SH1c.init, IH1c.init, DH1c.init, SH1n.init, IH1n.init, DH1n.init))
  popL2 <- sum(c(SL2c.init, IL2c.init, DL2c.init, SL2n.init, IL2n.init, DL2n.init))
  popH2 <- sum(c(SH2c.init, IH2c.init, DH2c.init, SH2n.init, IH2n.init, DH2n.init))
  
  
  # Balancing total discharge and admission rates
  rSH <- rSL/rH # need to do it this way because multiplier is meant to multiply the LOS not the rate!
  rIL <- rSL/rI
  rIH <- rSL/(rI*rH)
  rDL <- rD*rIL
  rDH <- rD*rIH
  hosp.discharge.rate <- (rSL*SL1c.init + rSH*SH1c.init + rIL*IL1c.init + rIH*IH1c.init + rDL*DL1c.init + rDH*DH1c.init +
                            rSL*SL1n.init + rSH*SH1n.init + rIL*IL1n.init + rIH*IH1n.init + rDL*DL1n.init + rDH*DH1n.init) 
  
  hospitalization.rate <- (hosp.discharge.rate)/(SL2c.init + rhoH*SH2c.init + rhoI*IL2c.init  + rhoH*rhoI*IH2c.init + rhoD*rhoI*DL2c.init + rhoD*rhoH*rhoI*DH2c.init +
                                                   SL2n.init + rhoH*SH2n.init + rhoI*IL2n.init  + rhoH*rhoI*IH2n.init + rhoD*rhoI*DL2n.init + rhoD*rhoH*rhoI*DH2n.init)
  aSL <- hospitalization.rate 
  aSH <- rhoH*hospitalization.rate
  aIL <- rhoI*hospitalization.rate
  aIH <- rhoH*rhoI*hospitalization.rate
  aDL <- rhoD*rhoI*hospitalization.rate
  aDH <- rhoD*rhoH*rhoI*hospitalization.rate
  
  # 8 different types of transmission
  beta1 <- old.beta1*(1 - int)
  
  L1infL1 <- beta1*(delta1 + (1 - delta1)*popL1/pop1)*(1/popL1)
  L2infL2 <- beta2*(delta2 + (1 - delta2)*popL2/pop2)*(1/popL2)
  L1infH1 <- beta1*((1 - delta1)*popL1/pop1)*(1/popL1)
  L2infH2 <- beta2*sigma*((1 - delta2)*popL2/pop2)*(1/popL2)
  H1infL1 <- beta1*((1 - delta1)*popH1/pop1)*(1/popH1)
  H2infL2 <- beta2*((1 - delta2)*popH2/pop2)*(1/popH2)
  H1infH1 <- beta1*(delta1 + (1 - delta1)*popH1/pop1)*(1/popH1)
  H2infH2 <- beta2*sigma*(delta2 + (1 - delta2)*popH2/pop2)*(1/popH2)
  
  
  IL1c_inf_IL1c <- SL1c.init*L1infL1
  IL1c_inf_IL1n <- SL1n.init*L1infL1
  IL1c_inf_IH1c <- SH1c.init*L1infH1
  IL1c_inf_IH1n <- SH1n.init*L1infH1
  IL1c_inf_IL2c <- 0
  IL1c_inf_IL2n <- 0
  IL1c_inf_IH2c <- 0
  IL1c_inf_IH2n <- 0
  IL1c_inf_DL1c <- 0
  IL1c_inf_DL1n <- 0
  IL1c_inf_DH1c <- 0
  IL1c_inf_DH1n <- 0
  IL1c_inf_DL2c <- 0
  IL1c_inf_DL2n <- 0
  IL1c_inf_DH2c <- 0
  IL1c_inf_DH2n <- 0
  
  IL1n_inf_IL1c <- SL1c.init*L1infL1
  IL1n_inf_IL1n <- SL1n.init*L1infL1
  IL1n_inf_IH1c <- SH1c.init*L1infH1
  IL1n_inf_IH1n <- SH1n.init*L1infH1
  IL1n_inf_IL2c <- 0
  IL1n_inf_IL2n <- 0
  IL1n_inf_IH2c <- 0
  IL1n_inf_IH2n <- 0
  IL1n_inf_DL1c <- 0
  IL1n_inf_DL1n <- 0
  IL1n_inf_DH1c <- 0
  IL1n_inf_DH1n <- 0
  IL1n_inf_DL2c <- 0
  IL1n_inf_DL2n <- 0
  IL1n_inf_DH2c <- 0
  IL1n_inf_DH2n <- 0
  
  IH1c_inf_IL1c <- SL1c.init*H1infL1
  IH1c_inf_IL1n <- SL1n.init*H1infL1
  IH1c_inf_IH1c <- SH1c.init*H1infH1
  IH1c_inf_IH1n <- SH1n.init*H1infH1
  IH1c_inf_IL2c <- 0
  IH1c_inf_IL2n <- 0
  IH1c_inf_IH2c <- 0
  IH1c_inf_IH2n <- 0
  IH1c_inf_DL1c <- 0
  IH1c_inf_DL1n <- 0
  IH1c_inf_DH1c <- 0
  IH1c_inf_DH1n <- 0
  IH1c_inf_DL2c <- 0
  IH1c_inf_DL2n <- 0
  IH1c_inf_DH2c <- 0
  IH1c_inf_DH2n <- 0
  
  IH1n_inf_IL1c <- SL1c.init*H1infL1
  IH1n_inf_IL1n <- SL1n.init*H1infL1
  IH1n_inf_IH1c <- SH1c.init*H1infH1
  IH1n_inf_IH1n <- SH1n.init*H1infH1
  IH1n_inf_IL2c <- 0
  IH1n_inf_IL2n <- 0
  IH1n_inf_IH2c <- 0
  IH1n_inf_IH2n <- 0
  IH1n_inf_DL1c <- 0
  IH1n_inf_DL1n <- 0
  IH1n_inf_DH1c <- 0
  IH1n_inf_DH1n <- 0
  IH1n_inf_DL2c <- 0
  IH1n_inf_DL2n <- 0
  IH1n_inf_DH2c <- 0
  IH1n_inf_DH2n <- 0
  
  IL2c_inf_IL1c <- 0
  IL2c_inf_IL1n <- 0
  IL2c_inf_IH1c <- 0
  IL2c_inf_IH1n <- 0
  IL2c_inf_IL2c <- SL2c.init*L2infL2
  IL2c_inf_IL2n <- SL2n.init*L2infL2
  IL2c_inf_IH2c <- SH2c.init*L2infH2
  IL2c_inf_IH2n <- SH2n.init*L2infH2
  IL2c_inf_DL1c <- 0
  IL2c_inf_DL1n <- 0
  IL2c_inf_DH1c <- 0
  IL2c_inf_DH1n <- 0
  IL2c_inf_DL2c <- 0
  IL2c_inf_DL2n <- 0
  IL2c_inf_DH2c <- 0
  IL2c_inf_DH2n <- 0
  
  IL2n_inf_IL1c <- 0
  IL2n_inf_IL1n <- 0
  IL2n_inf_IH1c <- 0
  IL2n_inf_IH1n <- 0
  IL2n_inf_IL2c <- SL2c.init*L2infL2
  IL2n_inf_IL2n <- SL2n.init*L2infL2
  IL2n_inf_IH2c <- SH2c.init*L2infH2
  IL2n_inf_IH2n <- SH2n.init*L2infH2
  IL2n_inf_DL1c <- 0
  IL2n_inf_DL1n <- 0
  IL2n_inf_DH1c <- 0
  IL2n_inf_DH1n <- 0
  IL2n_inf_DL2c <- 0
  IL2n_inf_DL2n <- 0
  IL2n_inf_DH2c <- 0
  IL2n_inf_DH2n <- 0
  
  IH2c_inf_IL1c <- 0
  IH2c_inf_IL1n <- 0
  IH2c_inf_IH1c <- 0
  IH2c_inf_IH1n <- 0
  IH2c_inf_IL2c <- SL2c.init*H2infL2
  IH2c_inf_IL2n <- SL2n.init*H2infL2
  IH2c_inf_IH2c <- SH2c.init*H2infH2
  IH2c_inf_IH2n <- SH2n.init*H2infH2
  IH2c_inf_DL1c <- 0
  IH2c_inf_DL1n <- 0
  IH2c_inf_DH1c <- 0
  IH2c_inf_DH1n <- 0
  IH2c_inf_DL2c <- 0
  IH2c_inf_DL2n <- 0
  IH2c_inf_DH2c <- 0
  IH2c_inf_DH2n <- 0
  
  IH2n_inf_IL1c <- 0
  IH2n_inf_IL1n <- 0
  IH2n_inf_IH1c <- 0
  IH2n_inf_IH1n <- 0
  IH2n_inf_IL2c <- SL2c.init*H2infL2
  IH2n_inf_IL2n <- SL2n.init*H2infL2
  IH2n_inf_IH2c <- SH2c.init*H2infH2
  IH2n_inf_IH2n <- SH2n.init*H2infH2
  IH2n_inf_DL1c <- 0
  IH2n_inf_DL1n <- 0
  IH2n_inf_DH1c <- 0
  IH2n_inf_DH1n <- 0
  IH2n_inf_DL2c <- 0
  IH2n_inf_DL2n <- 0
  IH2n_inf_DH2c <- 0
  IH2n_inf_DH2n <- 0
  
  
  DL1c_inf_IL1c <- SL1c.init*L1infL1
  DL1c_inf_IL1n <- SL1n.init*L1infL1
  DL1c_inf_IH1c <- SH1c.init*L1infH1
  DL1c_inf_IH1n <- SH1n.init*L1infH1
  DL1c_inf_IL2c <- 0
  DL1c_inf_IL2n <- 0
  DL1c_inf_IH2c <- 0
  DL1c_inf_IH2n <- 0
  DL1c_inf_DL1c <- 0
  DL1c_inf_DL1n <- 0
  DL1c_inf_DH1c <- 0
  DL1c_inf_DH1n <- 0
  DL1c_inf_DL2c <- 0
  DL1c_inf_DL2n <- 0
  DL1c_inf_DH2c <- 0
  DL1c_inf_DH2n <- 0
  
  DL1n_inf_IL1c <- SL1c.init*L1infL1
  DL1n_inf_IL1n <- SL1n.init*L1infL1
  DL1n_inf_IH1c <- SH1c.init*L1infH1
  DL1n_inf_IH1n <- SH1n.init*L1infH1
  DL1n_inf_IL2c <- 0
  DL1n_inf_IL2n <- 0
  DL1n_inf_IH2c <- 0
  DL1n_inf_IH2n <- 0
  DL1n_inf_DL1c <- 0
  DL1n_inf_DL1n <- 0
  DL1n_inf_DH1c <- 0
  DL1n_inf_DH1n <- 0
  DL1n_inf_DL2c <- 0
  DL1n_inf_DL2n <- 0
  DL1n_inf_DH2c <- 0
  DL1n_inf_DH2n <- 0
  
  DH1c_inf_IL1c <- SL1c.init*H1infL1
  DH1c_inf_IL1n <- SL1n.init*H1infL1
  DH1c_inf_IH1c <- SH1c.init*H1infH1
  DH1c_inf_IH1n <- SH1n.init*H1infH1
  DH1c_inf_IL2c <- 0
  DH1c_inf_IL2n <- 0
  DH1c_inf_IH2c <- 0
  DH1c_inf_IH2n <- 0
  DH1c_inf_DL1c <- 0
  DH1c_inf_DL1n <- 0
  DH1c_inf_DH1c <- 0
  DH1c_inf_DH1n <- 0
  DH1c_inf_DL2c <- 0
  DH1c_inf_DL2n <- 0
  DH1c_inf_DH2c <- 0
  DH1c_inf_DH2n <- 0
  
  DH1n_inf_IL1c <- SL1c.init*H1infL1
  DH1n_inf_IL1n <- SL1n.init*H1infL1
  DH1n_inf_IH1c <- SH1c.init*H1infH1
  DH1n_inf_IH1n <- SH1n.init*H1infH1
  DH1n_inf_IL2c <- 0
  DH1n_inf_IL2n <- 0
  DH1n_inf_IH2c <- 0
  DH1n_inf_IH2n <- 0
  DH1n_inf_DL1c <- 0
  DH1n_inf_DL1n <- 0
  DH1n_inf_DH1c <- 0
  DH1n_inf_DH1n <- 0
  DH1n_inf_DL2c <- 0
  DH1n_inf_DL2n <- 0
  DH1n_inf_DH2c <- 0
  DH1n_inf_DH2n <- 0
  
  DL2c_inf_IL1c <- 0
  DL2c_inf_IL1n <- 0
  DL2c_inf_IH1c <- 0
  DL2c_inf_IH1n <- 0
  DL2c_inf_IL2c <- SL2c.init*L2infL2
  DL2c_inf_IL2n <- SL2n.init*L2infL2
  DL2c_inf_IH2c <- SH2c.init*L2infH2
  DL2c_inf_IH2n <- SH2n.init*L2infH2
  DL2c_inf_DL1c <- 0
  DL2c_inf_DL1n <- 0
  DL2c_inf_DH1c <- 0
  DL2c_inf_DH1n <- 0
  DL2c_inf_DL2c <- 0
  DL2c_inf_DL2n <- 0
  DL2c_inf_DH2c <- 0
  DL2c_inf_DH2n <- 0
  
  DL2n_inf_IL1c <- 0
  DL2n_inf_IL1n <- 0
  DL2n_inf_IH1c <- 0
  DL2n_inf_IH1n <- 0
  DL2n_inf_IL2c <- SL2c.init*L2infL2
  DL2n_inf_IL2n <- SL2n.init*L2infL2
  DL2n_inf_IH2c <- SH2c.init*L2infH2
  DL2n_inf_IH2n <- SH2n.init*L2infH2
  DL2n_inf_DL1c <- 0
  DL2n_inf_DL1n <- 0
  DL2n_inf_DH1c <- 0
  DL2n_inf_DH1n <- 0
  DL2n_inf_DL2c <- 0
  DL2n_inf_DL2n <- 0
  DL2n_inf_DH2c <- 0
  DL2n_inf_DH2n <- 0
  
  DH2c_inf_IL1c <- 0
  DH2c_inf_IL1n <- 0
  DH2c_inf_IH1c <- 0
  DH2c_inf_IH1n <- 0
  DH2c_inf_IL2c <- SL2c.init*H2infL2
  DH2c_inf_IL2n <- SL2n.init*H2infL2
  DH2c_inf_IH2c <- SH2c.init*H2infH2
  DH2c_inf_IH2n <- SH2n.init*H2infH2
  DH2c_inf_DL1c <- 0
  DH2c_inf_DL1n <- 0
  DH2c_inf_DH1c <- 0
  DH2c_inf_DH1n <- 0
  DH2c_inf_DL2c <- 0
  DH2c_inf_DL2n <- 0
  DH2c_inf_DH2c <- 0
  DH2c_inf_DH2n <- 0
  
  DH2n_inf_IL1c <- 0
  DH2n_inf_IL1n <- 0
  DH2n_inf_IH1c <- 0
  DH2n_inf_IH1n <- 0
  DH2n_inf_IL2c <- SL2c.init*H2infL2
  DH2n_inf_IL2n <- SL2n.init*H2infL2
  DH2n_inf_IH2c <- SH2c.init*H2infH2
  DH2n_inf_IH2n <- SH2n.init*H2infH2
  DH2n_inf_DL1c <- 0
  DH2n_inf_DL1n <- 0
  DH2n_inf_DH1c <- 0
  DH2n_inf_DH1n <- 0
  DH2n_inf_DL2c <- 0
  DH2n_inf_DL2n <- 0
  DH2n_inf_DH2c <- 0
  DH2n_inf_DH2n <- 0
  
  Transmission.matrix = matrix(c(
    IL1c_inf_IL1c, IL1c_inf_IL1n, IL1c_inf_IH1c, IL1c_inf_IH1n, IL1c_inf_IL2c, IL1c_inf_IL2n, IL1c_inf_IH2c, IL1c_inf_IH2n, IL1c_inf_DL1c, IL1c_inf_DL1n, IL1c_inf_DH1c, IL1c_inf_DH1n, IL1c_inf_DL2c, IL1c_inf_DL2n, IL1c_inf_DH2c, IL1c_inf_DH2n,
    IL1n_inf_IL1c, IL1n_inf_IL1n, IL1n_inf_IH1c, IL1n_inf_IH1n, IL1n_inf_IL2c, IL1n_inf_IL2n, IL1n_inf_IH2c, IL1n_inf_IH2n, IL1n_inf_DL1c, IL1n_inf_DL1n, IL1n_inf_DH1c, IL1n_inf_DH1n, IL1n_inf_DL2c, IL1n_inf_DL2n, IL1n_inf_DH2c, IL1n_inf_DH2n,
    IH1c_inf_IL1c, IH1c_inf_IL1n, IH1c_inf_IH1c, IH1c_inf_IH1n, IH1c_inf_IL2c, IH1c_inf_IL2n, IH1c_inf_IH2c, IH1c_inf_IH2n, IH1c_inf_DL1c, IH1c_inf_DL1n, IH1c_inf_DH1c, IH1c_inf_DH1n, IH1c_inf_DL2c, IH1c_inf_DL2n, IH1c_inf_DH2c, IH1c_inf_DH2n,
    IH1n_inf_IL1c, IH1n_inf_IL1n, IH1n_inf_IH1c, IH1n_inf_IH1n, IH1n_inf_IL2c, IH1n_inf_IL2n, IH1n_inf_IH2c, IH1n_inf_IH2n, IH1n_inf_DL1c, IH1n_inf_DL1n, IH1n_inf_DH1c, IH1n_inf_DH1n, IH1n_inf_DL2c, IH1n_inf_DL2n, IH1n_inf_DH2c, IH1n_inf_DH2n,
    IL2c_inf_IL1c, IL2c_inf_IL1n, IL2c_inf_IH1c, IL2c_inf_IH1n, IL2c_inf_IL2c, IL2c_inf_IL2n, IL2c_inf_IH2c, IL2c_inf_IH2n, IL2c_inf_DL1c, IL2c_inf_DL1n, IL2c_inf_DH1c, IL2c_inf_DH1n, IL2c_inf_DL2c, IL2c_inf_DL2n, IL2c_inf_DH2c, IL2c_inf_DH2n,
    IL2n_inf_IL1c, IL2n_inf_IL1n, IL2n_inf_IH1c, IL2n_inf_IH1n, IL2n_inf_IL2c, IL2n_inf_IL2n, IL2n_inf_IH2c, IL2n_inf_IH2n, IL2n_inf_DL1c, IL2n_inf_DL1n, IL2n_inf_DH1c, IL2n_inf_DH1n, IL2n_inf_DL2c, IL2n_inf_DL2n, IL2n_inf_DH2c, IL2n_inf_DH2n,
    IH2c_inf_IL1c, IH2c_inf_IL1n, IH2c_inf_IH1c, IH2c_inf_IH1n, IH2c_inf_IL2c, IH2c_inf_IL2n, IH2c_inf_IH2c, IH2c_inf_IH2n, IH2c_inf_DL1c, IH2c_inf_DL1n, IH2c_inf_DH1c, IH2c_inf_DH1n, IH2c_inf_DL2c, IH2c_inf_DL2n, IH2c_inf_DH2c, IH2c_inf_DH2n,
    IH2c_inf_IL1c, IH2n_inf_IL1n, IH2n_inf_IH1c, IH2n_inf_IH1n, IH2n_inf_IL2c, IH2n_inf_IL2n, IH2n_inf_IH2c, IH2n_inf_IH2n, IH2n_inf_DL1c, IH2n_inf_DL1n, IH2n_inf_DH1c, IH2n_inf_DH1n, IH2n_inf_DL2c, IH2n_inf_DL2n, IH2n_inf_DH2c, IH2n_inf_DH2n,
    DL1c_inf_IL1c, DL1c_inf_IL1n, DL1c_inf_IH1c, DL1c_inf_IH1n, DL1c_inf_IL2c, DL1c_inf_IL2n, DL1c_inf_IH2c, DL1c_inf_IH2n, DL1c_inf_DL1c, DL1c_inf_DL1n, DL1c_inf_DH1c, DL1c_inf_DH1n, DL1c_inf_DL2c, DL1c_inf_DL2n, DL1c_inf_DH2c, DL1c_inf_DH2n,
    DL1n_inf_IL1c, DL1n_inf_IL1n, DL1n_inf_IH1c, DL1n_inf_IH1n, DL1n_inf_IL2c, DL1n_inf_IL2n, DL1n_inf_IH2c, DL1n_inf_IH2n, DL1n_inf_DL1c, DL1n_inf_DL1n, DL1n_inf_DH1c, DL1n_inf_DH1n, DL1n_inf_DL2c, DL1n_inf_DL2n, DL1n_inf_DH2c, DL1n_inf_DH2n,
    DH1c_inf_IL1c, DH1c_inf_IL1n, DH1c_inf_IH1c, DH1c_inf_IH1n, DH1c_inf_IL2c, DH1c_inf_IL2n, DH1c_inf_IH2c, DH1c_inf_IH2n, DH1c_inf_DL1c, DH1c_inf_DL1n, DH1c_inf_DH1c, DH1c_inf_DH1n, DH1c_inf_DL2c, DH1c_inf_DL2n, DH1c_inf_DH2c, DH1c_inf_DH2n,
    DH1n_inf_IL1c, DH1n_inf_IL1n, DH1n_inf_IH1c, DH1n_inf_IH1n, DH1n_inf_IL2c, DH1n_inf_IL2n, DH1n_inf_IH2c, DH1n_inf_IH2n, DH1n_inf_DL1c, DH1n_inf_DL1n, DH1n_inf_DH1c, DH1n_inf_DH1n, DH1n_inf_DL2c, DH1n_inf_DL2n, DH1n_inf_DH2c, DH1n_inf_DH2n,
    DL2c_inf_IL1c, DL2c_inf_IL1n, DL2c_inf_IH1c, DL2c_inf_IH1n, DL2c_inf_IL2c, DL2c_inf_IL2n, DL2c_inf_IH2c, DL2c_inf_IH2n, DL2c_inf_DL1c, DL2c_inf_DL1n, DL2c_inf_DH1c, DL2c_inf_DH1n, DL2c_inf_DL2c, DL2c_inf_DL2n, DL2c_inf_DH2c, DL2c_inf_DH2n,
    DL2n_inf_IL1c, DL2n_inf_IL1n, DL2n_inf_IH1c, DL2n_inf_IH1n, DL2n_inf_IL2c, DL2n_inf_IL2n, DL2n_inf_IH2c, DL2n_inf_IH2n, DL2n_inf_DL1c, DL2n_inf_DL1n, DL2n_inf_DH1c, DL2n_inf_DH1n, DL2n_inf_DL2c, DL2n_inf_DL2n, DL2n_inf_DH2c, DL2n_inf_DH2n,
    DH2c_inf_IL1c, DH2c_inf_IL1n, DH2c_inf_IH1c, DH2c_inf_IH1n, DH2c_inf_IL2c, DH2c_inf_IL2n, DH2c_inf_IH2c, DH2c_inf_IH2n, DH2c_inf_DL1c, DH2c_inf_DL1n, DH2c_inf_DH1c, DH2c_inf_DH1n, DH2c_inf_DL2c, DH2c_inf_DL2n, DH2c_inf_DH2c, DH2c_inf_DH2n,
    DH2n_inf_IL1c, DH2n_inf_IL1n, DH2n_inf_IH1c, DH2n_inf_IH1n, DH2n_inf_IL2c, DH2n_inf_IL2n, DH2n_inf_IH2c, DH2n_inf_IH2n, DH2n_inf_DL1c, DH2n_inf_DL1n, DH2n_inf_DH1c, DH2n_inf_DH1n, DH2n_inf_DL2c, DH2n_inf_DL2n, DH2n_inf_DH2c, DH2n_inf_DH2n),
    nrow = 16, byrow = FALSE)
  
  IL1c_to_IL1c <- - (gamma1c + rIL + alphaIL1)
  IL1c_to_IL1n <- 0
  IL1c_to_IH1c <- 0
  IL1c_to_IH1n <- 0
  IL1c_to_IL2c <- rIL
  IL1c_to_IL2n <- 0
  IL1c_to_IH2c <- 0
  IL1c_to_IH2n <- 0
  IL1c_to_DL1c <- alphaIL1
  IL1c_to_DL1n <- 0
  IL1c_to_DH1c <- 0
  IL1c_to_DH1n <- 0
  IL1c_to_DL2c <- 0
  IL1c_to_DL2n <- 0
  IL1c_to_DH2c <- 0
  IL1c_to_DH2n <- 0
  
  IL1n_to_IL1c <- 0
  IL1n_to_IL1n <- - (gamma1n + rIL + alphaIL1)
  IL1n_to_IH1c <- 0
  IL1n_to_IH1n <- 0
  IL1n_to_IL2c <- 0
  IL1n_to_IL2n <- rIL
  IL1n_to_IH2c <- 0
  IL1n_to_IH2n <- 0
  IL1n_to_DL1c <- 0
  IL1n_to_DL1n <- alphaIL1
  IL1n_to_DH1c <- 0
  IL1n_to_DH1n <- 0
  IL1n_to_DL2c <- 0
  IL1n_to_DL2n <- 0
  IL1n_to_DH2c <- 0
  IL1n_to_DH2n <- 0
  
  IH1c_to_IL1c <- 0
  IH1c_to_IL1n <- 0
  IH1c_to_IH1c <- - (gamma1c + rIH + alphaIH1)
  IH1c_to_IH1n <- 0
  IH1c_to_IL2c <- 0
  IH1c_to_IL2n <- 0
  IH1c_to_IH2c <- rIH
  IH1c_to_IH2n <- 0
  IH1c_to_DL1c <- 0
  IH1c_to_DL1n <- 0
  IH1c_to_DH1c <- alphaIH1
  IH1c_to_DH1n <- 0
  IH1c_to_DL2c <- 0
  IH1c_to_DL2n <- 0
  IH1c_to_DH2c <- 0
  IH1c_to_DH2n <- 0
  
  IH1n_to_IL1c <- 0
  IH1n_to_IL1n <- 0
  IH1n_to_IH1c <- 0
  IH1n_to_IH1n <- - (gamma1n + rIH + alphaIH1)
  IH1n_to_IL2c <- 0
  IH1n_to_IL2n <- 0
  IH1n_to_IH2c <- 0
  IH1n_to_IH2n <- rIH
  IH1n_to_DL1c <- 0
  IH1n_to_DL1n <- 0
  IH1n_to_DH1c <- 0
  IH1n_to_DH1n <- alphaIH1
  IH1n_to_DL2c <- 0
  IH1n_to_DL2n <- 0
  IH1n_to_DH2c <- 0
  IH1n_to_DH2n <- 0
  
  IL2c_to_IL1c <- aIL
  IL2c_to_IL1n <- 0
  IL2c_to_IH1c <- 0
  IL2c_to_IH1n <- 0
  IL2c_to_IL2c <- - (gamma2c + aIL + alphaIL2)
  IL2c_to_IL2n <- 0
  IL2c_to_IH2c <- 0
  IL2c_to_IH2n <- 0
  IL2c_to_DL1c <- 0
  IL2c_to_DL1n <- 0
  IL2c_to_DH1c <- 0
  IL2c_to_DH1n <- 0
  IL2c_to_DL2c <- alphaIL2
  IL2c_to_DL2n <- 0
  IL2c_to_DH2c <- 0
  IL2c_to_DH2n <- 0
  
  IL2n_to_IL1c <- 0
  IL2n_to_IL1n <- aIL
  IL2n_to_IH1c <- 0
  IL2n_to_IH1n <- 0
  IL2n_to_IL2c <- 0
  IL2n_to_IL2n <- - (gamma2n + aIL + alphaIL2)
  IL2n_to_IH2c <- 0
  IL2n_to_IH2n <- 0
  IL2n_to_DL1c <- 0
  IL2n_to_DL1n <- 0
  IL2n_to_DH1c <- 0
  IL2n_to_DH1n <- 0
  IL2n_to_DL2c <- 0
  IL2n_to_DL2n <- alphaIL2
  IL2n_to_DH2c <- 0
  IL2n_to_DH2n <- 0
  
  IH2c_to_IL1c <- 0
  IH2c_to_IL1n <- 0
  IH2c_to_IH1c <- aIH
  IH2c_to_IH1n <- 0
  IH2c_to_IL2c <- 0
  IH2c_to_IL2n <- 0
  IH2c_to_IH2c <- - (gamma2c + aIH + alphaIH2)
  IH2c_to_IH2n <- 0
  IH2c_to_DL1c <- 0
  IH2c_to_DL1n <- 0
  IH2c_to_DH1c <- 0
  IH2c_to_DH1n <- 0
  IH2c_to_DL2c <- 0
  IH2c_to_DL2n <- 0
  IH2c_to_DH2c <- alphaIH2
  IH2c_to_DH2n <- 0
  
  IH2n_to_IL1c <- 0
  IH2n_to_IL1n <- 0
  IH2n_to_IH1c <- 0
  IH2n_to_IH1n <- aIH
  IH2n_to_IL2c <- 0
  IH2n_to_IL2n <- 0
  IH2n_to_IH2c <- 0
  IH2n_to_IH2n <- - (gamma2n + aIH + alphaIH2)
  IH2n_to_DL1c <- 0
  IH2n_to_DL1n <- 0
  IH2n_to_DH1c <- 0
  IH2n_to_DH1n <- 0
  IH2n_to_DL2c <- 0
  IH2n_to_DL2n <- 0
  IH2n_to_DH2c <- 0
  IH2n_to_DH2n <- alphaIH2
  
  DL1c_to_IL1c <- gammaD
  DL1c_to_IL1n <- 0
  DL1c_to_IH1c <- 0
  DL1c_to_IH1n <- 0
  DL1c_to_IL2c <- 0
  DL1c_to_IL2n <- 0
  DL1c_to_IH2c <- 0
  DL1c_to_IH2n <- 0
  DL1c_to_DL1c <- - (gammaD + rDL)
  DL1c_to_DL1n <- 0
  DL1c_to_DH1c <- 0
  DL1c_to_DH1n <- 0
  DL1c_to_DL2c <- rDL
  DL1c_to_DL2n <- 0
  DL1c_to_DH2c <- 0
  DL1c_to_DH2n <- 0
  
  DL1n_to_IL1c <- 0
  DL1n_to_IL1n <- gammaD
  DL1n_to_IH1c <- 0
  DL1n_to_IH1n <- 0
  DL1n_to_IL2c <- 0
  DL1n_to_IL2n <- 0
  DL1n_to_IH2c <- 0
  DL1n_to_IH2n <- 0
  DL1n_to_DL1c <- 0
  DL1n_to_DL1n <- - (gammaD + rDL)
  DL1n_to_DH1c <- 0
  DL1n_to_DH1n <- 0
  DL1n_to_DL2c <- 0
  DL1n_to_DL2n <- rDL
  DL1n_to_DH2c <- 0
  DL1n_to_DH2n <- 0
  
  DH1c_to_IL1c <- 0
  DH1c_to_IL1n <- 0
  DH1c_to_IH1c <- gammaD
  DH1c_to_IH1n <- 0
  DH1c_to_IL2c <- 0
  DH1c_to_IL2n <- 0
  DH1c_to_IH2c <- 0
  DH1c_to_IH2n <- 0
  DH1c_to_DL1c <- 0
  DH1c_to_DL1n <- 0
  DH1c_to_DH1c <- - (gammaD + rDH)
  DH1c_to_DH1n <- 0
  DH1c_to_DL2c <- 0
  DH1c_to_DL2n <- 0
  DH1c_to_DH2c <- rDH
  DH1c_to_DH2n <- 0
  
  DH1n_to_IL1c <- 0
  DH1n_to_IL1n <- 0
  DH1n_to_IH1c <- 0
  DH1n_to_IH1n <- gammaD
  DH1n_to_IL2c <- 0
  DH1n_to_IL2n <- 0
  DH1n_to_IH2c <- 0
  DH1n_to_IH2n <- 0
  DH1n_to_DL1c <- 0
  DH1n_to_DL1n <- 0
  DH1n_to_DH1c <- 0
  DH1n_to_DH1n <- - (gammaD + rDH)
  DH1n_to_DL2c <- 0
  DH1n_to_DL2n <- 0
  DH1n_to_DH2c <- 0
  DH1n_to_DH2n <- rDH
  
  DL2c_to_IL1c <- 0
  DL2c_to_IL1n <- 0
  DL2c_to_IH1c <- 0
  DL2c_to_IH1n <- 0
  DL2c_to_IL2c <- gammaD
  DL2c_to_IL2n <- 0
  DL2c_to_IH2c <- 0
  DL2c_to_IH2n <- 0
  DL2c_to_DL1c <- aDL
  DL2c_to_DL1n <- 0
  DL2c_to_DH1c <- 0
  DL2c_to_DH1n <- 0
  DL2c_to_DL2c <- - (gammaD + aDL) 
  DL2c_to_DL2n <- 0
  DL2c_to_DH2c <- 0
  DL2c_to_DH2n <- 0
  
  DL2n_to_IL1c <- 0
  DL2n_to_IL1n <- 0
  DL2n_to_IH1c <- 0
  DL2n_to_IH1n <- 0
  DL2n_to_IL2c <- 0
  DL2n_to_IL2n <- gammaD
  DL2n_to_IH2c <- 0
  DL2n_to_IH2n <- 0
  DL2n_to_DL1c <- 0
  DL2n_to_DL1n <- aDL
  DL2n_to_DH1c <- 0
  DL2n_to_DH1n <- 0
  DL2n_to_DL2c <- 0
  DL2n_to_DL2n <- - (gammaD + aDL) 
  DL2n_to_DH2c <- 0
  DL2n_to_DH2n <- 0
  
  DH2c_to_IL1c <- 0
  DH2c_to_IL1n <- 0
  DH2c_to_IH1c <- 0
  DH2c_to_IH1n <- 0
  DH2c_to_IL2c <- 0
  DH2c_to_IL2n <- 0
  DH2c_to_IH2c <- gammaD
  DH2c_to_IH2n <- 0
  DH2c_to_DL1c <- 0
  DH2c_to_DL1n <- 0
  DH2c_to_DH1c <- aDH
  DH2c_to_DH1n <- 0
  DH2c_to_DL2c <- 0
  DH2c_to_DL2n <- 0
  DH2c_to_DH2c <- - (gammaD + aDH) 
  DH2c_to_DH2n <- 0
  
  DH2n_to_IL1c <- 0
  DH2n_to_IL1n <- 0
  DH2n_to_IH1c <- 0
  DH2n_to_IH1n <- 0
  DH2n_to_IL2c <- 0
  DH2n_to_IL2n <- 0
  DH2n_to_IH2c <- 0
  DH2n_to_IH2n <- gammaD
  DH2n_to_DL1c <- 0
  DH2n_to_DL1n <- 0
  DH2n_to_DH1c <- 0
  DH2n_to_DH1n <- aDH
  DH2n_to_DL2c <- 0
  DH2n_to_DL2n <- 0
  DH2n_to_DH2c <- 0
  DH2n_to_DH2n <- - (gammaD + aDH) 
  
  Transition.matrix = matrix(c(
    IL1c_to_IL1c, IL1c_to_IL1n, IL1c_to_IH1c, IL1c_to_IH1n, IL1c_to_IL2c, IL1c_to_IL2n, IL1c_to_IH2c, IL1c_to_IH2n, IL1c_to_DL1c, IL1c_to_DL1n, IL1c_to_DH1c, IL1c_to_DH1n, IL1c_to_DL2c, IL1c_to_DL2n, IL1c_to_DH2c, IL1c_to_DH2n,
    IL1n_to_IL1c, IL1n_to_IL1n, IL1n_to_IH1c, IL1n_to_IH1n, IL1n_to_IL2c, IL1n_to_IL2n, IL1n_to_IH2c, IL1n_to_IH2n, IL1n_to_DL1c, IL1n_to_DL1n, IL1n_to_DH1c, IL1n_to_DH1n, IL1n_to_DL2c, IL1n_to_DL2n, IL1n_to_DH2c, IL1n_to_DH2n,
    IH1c_to_IL1c, IH1c_to_IL1n, IH1c_to_IH1c, IH1c_to_IH1n, IH1c_to_IL2c, IH1c_to_IL2n, IH1c_to_IH2c, IH1c_to_IH2n, IH1c_to_DL1c, IH1c_to_DL1n, IH1c_to_DH1c, IH1c_to_DH1n, IH1c_to_DL2c, IH1c_to_DL2n, IH1c_to_DH2c, IH1c_to_DH2n,
    IH1n_to_IL1c, IH1n_to_IL1n, IH1n_to_IH1c, IH1n_to_IH1n, IH1n_to_IL2c, IH1n_to_IL2n, IH1n_to_IH2c, IH1n_to_IH2n, IH1n_to_DL1c, IH1n_to_DL1n, IH1n_to_DH1c, IH1n_to_DH1n, IH1n_to_DL2c, IH1n_to_DL2n, IH1n_to_DH2c, IH1n_to_DH2n,
    IL2c_to_IL1c, IL2c_to_IL1n, IL2c_to_IH1c, IL2c_to_IH1n, IL2c_to_IL2c, IL2c_to_IL2n, IL2c_to_IH2c, IL2c_to_IH2n, IL2c_to_DL1c, IL2c_to_DL1n, IL2c_to_DH1c, IL2c_to_DH1n, IL2c_to_DL2c, IL2c_to_DL2n, IL2c_to_DH2c, IL2c_to_DH2n,
    IL2n_to_IL1c, IL2n_to_IL1n, IL2n_to_IH1c, IL2n_to_IH1n, IL2n_to_IL2c, IL2n_to_IL2n, IL2n_to_IH2c, IL2n_to_IH2n, IL2n_to_DL1c, IL2n_to_DL1n, IL2n_to_DH1c, IL2n_to_DH1n, IL2n_to_DL2c, IL2n_to_DL2n, IL2n_to_DH2c, IL2n_to_DH2n,
    IH2c_to_IL1c, IH2c_to_IL1n, IH2c_to_IH1c, IH2c_to_IH1n, IH2c_to_IL2c, IH2c_to_IL2n, IH2c_to_IH2c, IH2c_to_IH2n, IH2c_to_DL1c, IH2c_to_DL1n, IH2c_to_DH1c, IH2c_to_DH1n, IH2c_to_DL2c, IH2c_to_DL2n, IH2c_to_DH2c, IH2c_to_DH2n,
    IH2n_to_IL1c, IH2n_to_IL1n, IH2n_to_IH1c, IH2n_to_IH1n, IH2n_to_IL2c, IH2n_to_IL2n, IH2n_to_IH2c, IH2n_to_IH2n, IH2n_to_DL1c, IH2n_to_DL1n, IH2n_to_DH1c, IH2n_to_DH1n, IH2n_to_DL2c, IH2n_to_DL2n, IH2n_to_DH2c, IH2n_to_DH2n,
    DL1c_to_IL1c, DL1c_to_IL1n, DL1c_to_IH1c, DL1c_to_IH1n, DL1c_to_IL2c, DL1c_to_IL2n, DL1c_to_IH2c, DL1c_to_IH2n, DL1c_to_DL1c, DL1c_to_DL1n, DL1c_to_DH1c, DL1c_to_DH1n, DL1c_to_DL2c, DL1c_to_DL2n, DL1c_to_DH2c, DL1c_to_DH2n,
    DL1n_to_IL1c, DL1n_to_IL1n, DL1n_to_IH1c, DL1n_to_IH1n, DL1n_to_IL2c, DL1n_to_IL2n, DL1n_to_IH2c, DL1n_to_IH2n, DL1n_to_DL1c, DL1n_to_DL1n, DL1n_to_DH1c, DL1n_to_DH1n, DL1n_to_DL2c, DL1n_to_DL2n, DL1n_to_DH2c, DL1n_to_DH2n,
    DH1c_to_IL1c, DH1c_to_IL1n, DH1c_to_IH1c, DH1c_to_IH1n, DH1c_to_IL2c, DH1c_to_IL2n, DH1c_to_IH2c, DH1c_to_IH2n, DH1c_to_DL1c, DH1c_to_DL1n, DH1c_to_DH1c, DH1c_to_DH1n, DH1c_to_DL2c, DH1c_to_DL2n, DH1c_to_DH2c, DH1c_to_DH2n,
    DH1n_to_IL1c, DH1n_to_IL1n, DH1n_to_IH1c, DH1n_to_IH1n, DH1n_to_IL2c, DH1n_to_IL2n, DH1n_to_IH2c, DH1n_to_IH2n, DH1n_to_DL1c, DH1n_to_DL1n, DH1n_to_DH1c, DH1n_to_DH1n, DH1n_to_DL2c, DH1n_to_DL2n, DH1n_to_DH2c, DH1n_to_DH2n,
    DL2c_to_IL1c, DL2c_to_IL1n, DL2c_to_IH1c, DL2c_to_IH1n, DL2c_to_IL2c, DL2c_to_IL2n, DL2c_to_IH2c, DL2c_to_IH2n, DL2c_to_DL1c, DL2c_to_DL1n, DL2c_to_DH1c, DL2c_to_DH1n, DL2c_to_DL2c, DL2c_to_DL2n, DL2c_to_DH2c, DL2c_to_DH2n,
    DL2n_to_IL1c, DL2n_to_IL1n, DL2n_to_IH1c, DL2n_to_IH1n, DL2n_to_IL2c, DL2n_to_IL2n, DL2n_to_IH2c, DL2n_to_IH2n, DL2n_to_DL1c, DL2n_to_DL1n, DL2n_to_DH1c, DL2n_to_DH1n, DL2n_to_DL2c, DL2n_to_DL2n, DL2n_to_DH2c, DL2n_to_DH2n,
    DH2c_to_IL1c, DH2c_to_IL1n, DH2c_to_IH1c, DH2c_to_IH1n, DH2c_to_IL2c, DH2c_to_IL2n, DH2c_to_IH2c, DH2c_to_IH2n, DH2c_to_DL1c, DH2c_to_DL1n, DH2c_to_DH1c, DH2c_to_DH1n, DH2c_to_DL2c, DH2c_to_DL2n, DH2c_to_DH2c, DH2c_to_DH2n,
    DH2n_to_IL1c, DH2n_to_IL1n, DH2n_to_IH1c, DH2n_to_IH1n, DH2n_to_IL2c, DH2n_to_IL2n, DH2n_to_IH2c, DH2n_to_IH2n, DH2n_to_DL1c, DH2n_to_DL1n, DH2n_to_DH1c, DH2n_to_DH1n, DH2n_to_DL2c, DH2n_to_DL2n, DH2n_to_DH2c, DH2n_to_DH2n),
    nrow = 16, byrow = FALSE)
  
  ngm <- (-Transmission.matrix) %*% (solve(Transition.matrix)) #solve for inverse, multiply by negative transmission matrix
  
  R0 <- max(as.double(eigen(ngm)$values))
  
  
  
  #R_hosp 
  #transmission matrix only includes hospital transmissions
  #transition matrix is the same as above
  
  Transmission.matrix.hosp = matrix(c(
    IL1c_inf_IL1c, IL1c_inf_IL1n, IL1c_inf_IH1c, IL1c_inf_IH1n, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    IL1n_inf_IL1c, IL1n_inf_IL1n, IL1n_inf_IH1c, IL1n_inf_IH1n, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    IH1c_inf_IL1c, IH1c_inf_IL1n, IH1c_inf_IH1c, IH1c_inf_IH1n, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    IH1n_inf_IL1c, IH1n_inf_IL1n, IH1n_inf_IH1c, IH1n_inf_IH1n, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    DL1c_inf_IL1c, DL1c_inf_IL1n, DL1c_inf_IH1c, DL1c_inf_IH1n, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    DL1n_inf_IL1c, DL1n_inf_IL1n, DL1n_inf_IH1c, DL1n_inf_IH1n, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    DH1c_inf_IL1c, DH1c_inf_IL1n, DH1c_inf_IH1c, DH1c_inf_IH1n, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    DH1n_inf_IL1c, DH1n_inf_IL1n, DH1n_inf_IH1c, DH1n_inf_IH1n, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    nrow = 16, byrow = FALSE)
  
  ngm_hosp <- (-Transmission.matrix.hosp) %*% (solve(Transition.matrix)) 
  
  R_hosp <- max(as.double(eigen(ngm_hosp)$values))
  
  #R_comm
  #transmission matrix only includes community transmissions
  #transition matrix is the same as above
  
  Transmission.matrix.comm = matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, IL2c_inf_IL2c, IL2c_inf_IL2n, IL2c_inf_IH2c, IL2c_inf_IH2n, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, IL2n_inf_IL2c, IL2n_inf_IL2n, IL2n_inf_IH2c, IL2n_inf_IH2n, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, IH2c_inf_IL2c, IH2c_inf_IL2n, IH2c_inf_IH2c, IH2c_inf_IH2n, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, IH2n_inf_IL2c, IH2n_inf_IL2n, IH2n_inf_IH2c, IH2n_inf_IH2n, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, DL2c_inf_IL2c, DL2c_inf_IL2n, DL2c_inf_IH2c, DL2c_inf_IH2n, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, DL2n_inf_IL2c, DL2n_inf_IL2n, DL2n_inf_IH2c, DL2n_inf_IH2n, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, DH2c_inf_IL2c, DH2c_inf_IL2n, DH2c_inf_IH2c, DH2c_inf_IH2n, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, DH2n_inf_IL2c, DH2n_inf_IL2n, DH2n_inf_IH2c, DH2n_inf_IH2n, 0, 0, 0, 0, 0, 0, 0, 0), 
    nrow = 16, byrow = FALSE)
  
  ngm_comm <- (-Transmission.matrix.comm) %*% (solve(Transition.matrix)) 
  
  R_comm <- max(as.double(eigen(ngm_comm)$values))
  
  
  list(R0 = R0, R_hosp = R_hosp, R_comm = R_comm, ngm = ngm)
  
}

# SSE FITTING BETAS ===============================================================
sse.rs.prev <- function(params0) {
  hosp.to.fit <- 0.034
  comm.to.fit <- 0.015
  
  old.beta1 <- exp(params0[1])
  beta2 <- exp(params0[2])
  
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
  no.int.equil <- runsteady(y = mod.init, func = modsim.dyn, parms = modpars, 
                            times = c(0, 10000000))
  
  eqprev.hosp.noint <- sum(no.int.equil$y[c(3, 4, 5, 6, 15, 16, 17, 18)])/N1.t0
  eqprev.comm.noint <- sum(no.int.equil$y[c(9, 10, 11, 12, 21, 22, 23, 24)])/N2.t0
  
  sse <- sum((eqprev.hosp.noint - hosp.to.fit)^2, (eqprev.comm.noint - comm.to.fit)^2)
}