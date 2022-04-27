# heatmaps for simulataneous interventions

source("model_functions.R")
require(deSolve)
require(rootSolve)
require(tidyverse)
require(sensitivity)
require(lhs)
require(cowplot)
require(viridis)

assumptions_df <- read.csv("assumptions_beta_fiting.csv")
assumptions_df <- assumptions_df[-2,] # remove heterogeneous assumption

beta.int.vec <- seq(0, 1, length.out = 5)
decol.int.vec <- seq(1/365, 1/5, length.out = 5) # baseline carriage length = 275 days

both.int.vec <- cbind(beta.int.vec, decol.int.vec)
colnames(both.int.vec) <- NULL

both.data <- as.data.frame(matrix(rep(NA, 6*150), nrow = 150))
names(both.data) <- names(both.data) <-c("beta.int", "decol.rate", "hospital", 
                                         "community", "cases.averted", "assumption")

j <- 1

for (i in 1:nrow(assumptions_df)) {
  
  sigma <- assumptions_df[i, "sigma"] 
  delta1 <- assumptions_df[i, "delta1"]
  delta2 <- assumptions_df[i, "delta2"]
  rSL <- assumptions_df[i, "rSL"] 
  rSH <- assumptions_df[i, "rSH"] 
  rIL <- assumptions_df[i, "rIL"] 
  rIH <- assumptions_df[i, "rIH"] 
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
  
  for (f in 1:length(beta.int.vec)) {
    
    new.beta.int <- beta.int.vec[f]
    
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
      
      # 2 No intervention--deSolve to find 5 year incidence
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
                       rSH, rIL, rIH, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = new.beta.int, 
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
      
      # 4 Intervention--runsteady to find equilibrium 
      modpars.int <- c(old.beta1, delta1, beta2, delta2, sigma, gamma1n, gamma1c, gamma2n, 
                       gamma2c, gammaD, alphaIL1, alphaIH1, alphaIL2, alphaIH2, rSL, 
                       rSH, rIL, rIH, rD, rhoH, rhoI, rhoD, int, time.int = 0, beta.int = new.beta.int, 
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
      
      both.data[j, 1] <- new.beta.int
      both.data[j, 2] <- decol.int
      both.data[j, 3] <- eqprev.hosp.int
      both.data[j, 4] <- eqprev.comm.int
      both.data[j, 5] <- cases.averted
      both.data[j, 6] <- paste(assumptions_df[i, "assumption"])
      
      j <- j + 1
      print(j)
   }
  }
}  
  
write.csv(both.data, "simul_int_data.csv", row.names = FALSE)

# HEATMAP OF CASES AVERTED  =========================================================
both.data <- read.csv("simul_int_data.csv")

both.data <- both.data %>% filter(assumption != "heterogen")
both.data$assumption <- gsub("homogen", "Homogeneous", both.data$assumption)
both.data$assumption <- gsub("baseline", "Age Group", both.data$assumption)
both.data$assumption <- gsub("Risk Group", "Age Group", both.data$assumption)
both.data$assumption <- gsub("inf_adm", "Colonized Admission", both.data$assumption)
both.data$assumption <- gsub("inf_los", "Colonized LOS", both.data$assumption)
both.data$assumption <- gsub("inf_both", "Colonized LOS + Admission", both.data$assumption)
both.data$assumption <- as.factor(both.data$assumption)

# smoothed with linear interpolation
p.1a <- ggplot(data = both.data %>% filter(assumption != "Homogeneous"), 
             aes(x = beta.int*100, y = decol.rate, color = cases.averted)) +
  geom_raster(aes(fill = cases.averted), interpolate = TRUE) +
  facet_wrap(~assumption, ncol = 2, scales = "free") +
  scale_fill_viridis_c(guide_legend(title = "Cases Averted")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100"),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     labels = c("0", "0.05", "0.1", "0.15", "0.20"),
                     expand = c(0,0)) +
  xlab(expression(paste("Transmission intervention efficacy % ", "(",phi,")"))) +
  ylab(expression(paste("Decolonization rate", " (", gamma[int],")"))) +
  theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank())
p.1a

# tiled w/o interpolation
p.2a <- ggplot(data = both.data %>% filter(assumption != "Homogeneous"), 
              aes(x = beta.int*100, y = decol.rate, fill = cases.averted)) +
  geom_tile(aes(fill = cases.averted)) +
  facet_wrap(~assumption, ncol = 2, scales = "free") +
  scale_fill_viridis_c(guide_legend(title = "Cases Averted")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100"),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     labels = c("0", "0.05", "0.1", "0.15", "0.20"),
                     expand = c(0,0)) +
  xlab(expression(paste("Transmission intervention efficacy % ", "(",phi,")"))) +
  ylab(expression(paste("Decolonization rate", " (", gamma[int],")"))) +
  theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank())
p.2a

# add contour lines
p.3a <- ggplot(data = both.data %>% filter(assumption != "Homogeneous"), 
              aes(x = beta.int*100, y = decol.rate, z = cases.averted)) +
  geom_tile(aes(fill = cases.averted)) +
  #stat_contour(bins = 10, color = "white") +
  facet_wrap(~assumption, ncol = 2, scales = "free") +
  scale_fill_viridis_c(guide_legend(title = "Cases Averted")) +
  scale_colour_gradientn(colours = c("grey","grey"),
                         values = c(0,100)) +
  # stat_contour(bins = 10, color = "white") +
  stat_contour(aes(colour = ..level..)) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100"),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     labels = c("0", "0.05", "0.1", "0.15", "0.20"),
                     expand = c(0,0)) +
  xlab(expression(paste("Transmission intervention efficacy % ", "(",phi,")"))) +
  ylab(expression(paste("Decolonization rate", " (", gamma[int],")"))) +
  theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank())
p.3a

p.3alab <- direct.label(p.3a, list("top.pieces",
                                            hjust = 0, 
                                            vjust = 0,
                                            cex = 0.75, color = "white"))
p.3alab


# add homogeneous example
p.4a <- ggplot(data = both.data, 
              aes(x = beta.int*100, y = decol.rate, z = cases.averted)) +
  geom_tile(aes(fill = cases.averted)) +
  stat_contour(bins = 10, color = "white") +
  facet_wrap(~assumption, ncol = 3, scales = "free") +
  scale_fill_viridis_c(guide_legend(title = "Cases Averted")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100"),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     labels = c("0", "0.05", "0.1", "0.15", "0.20"),
                     expand = c(0,0)) +
  xlab(expression(paste("Transmission intervention efficacy % ", "(",phi,")"))) +
  ylab(expression(paste("Decolonization rate", " (", gamma[int],")"))) +
  theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank()) 
p.4a

pdf(file = "simul_int_cases_heatmap.pdf", width = 8, height = 5.5, useDingbats=F)
p.1a
p.2a
p.3a
p.4a
dev.off()


# add contour lines
p.3basic <- ggplot(data = both.data, 
               aes(x = beta.int*100, y = decol.rate, z = cases.averted/8223*100)) +
  geom_tile(aes(fill = cases.averted/8223*100)) +
  facet_wrap(~assumption, ncol = 2, scales = "free") +
  scale_fill_viridis_c(guide_legend(title = "% Cases Averted")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100"),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     labels = c("0", "0.05", "0.1", "0.15", "0.20"),
                     expand = c(0,0)) +
  xlab(expression(paste("Transmission intervention efficacy % ", "(",phi,")"))) +
  ylab(expression(paste("Decolonization rate", " (", gamma[int],")"))) +
  theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank())
p.3basic

p.3both <- ggplot(data = both.data, 
                   aes(x = beta.int*100, y = decol.rate, z = cases.averted/8223*100)) +
  geom_tile(aes(fill = cases.averted/8223*100)) +
  #stat_contour(bins = 10, color = "white") +
  facet_wrap(~assumption, ncol = 2, scales = "free") +
  scale_fill_viridis_c(guide_legend(title = "% Cases Averted")) +
  scale_colour_gradientn(colours = c("grey","grey"),
                         values = c(0,100)) +
  # stat_contour(bins = 10, color = "white") +
  stat_contour(aes(colour = ..level..)) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100"),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     labels = c("0", "0.05", "0.1", "0.15", "0.20"),
                     expand = c(0,0)) +
  xlab(expression(paste("Transmission intervention efficacy % ", "(",phi,")"))) +
  ylab(expression(paste("Decolonization rate", " (", gamma[int],")"))) +
  theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank())
p.3both


p3leg <- get_legend(p.3basic)

p.3.withlabs <- direct.label(p.3both, list("top.pieces",
                                   hjust = 0, 
                                   vjust = 0,
                                   cex = 0.75, color = "white"))


fullheatmap <- plot_grid(p.3.withlabs, p3leg, ncol = 2,
                rel_widths = c(1.5, .31),
                rel_heights = c(1.75, .31))
fullheatmap


pdf(file = "simul_int_cases_heatmap.pdf", width = 8.1, height = 8, useDingbats=F)
fullheatmap
dev.off()

# save eps for submission
save_plot("FigS3.eps", fullheatmap,
          base_height = 8, base_width = 8.1)



# HEATMAP OF EQUILIBRIUM PREVALENCES  =========================================================

# hospital equilibrium prevalence with intervention
p.hosp <- ggplot(data = both.data %>% filter(assumption != "Homogeneous"), 
                 aes(x = beta.int*100, y = decol.rate, z = hospital*100)) +
  geom_raster(aes(fill = hospital*100)) +
  #stat_contour(bins = 10, color = "white") +
  facet_wrap(~assumption, ncol = 2, scales = "free") +
  scale_fill_viridis(option = "A", guide_legend(title = "Hospital prevalence \nwith intervention")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100"),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     labels = c("0", "0.05", "0.1", "0.15", "0.20"),
                     expand = c(0,0)) +
  xlab(expression(paste("Transmission intervention efficacy % ", "(",phi,")"))) +
  ylab(expression(paste("Decolonization rate", " (", gamma[int],")"))) +
  theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank()) +
  ggtitle("Equilibrium prevalence in the hospital \nafter an intervention")
p.hosp

# community equilibrium prevalence with intervention
p.comm <- ggplot(data = both.data %>% filter(assumption != "Homogeneous"), 
                 aes(x = beta.int*100, y = decol.rate, z = community*100)) +
  geom_raster(aes(fill = community*100)) +
  #stat_contour(bins = 10, color = "white") +
  facet_wrap(~assumption, ncol = 2, scales = "free") +
  scale_fill_viridis(option = "A", guide_legend(title = "Community prevalence \nwith intervention")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100"),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     labels = c("0", "0.05", "0.1", "0.15", "0.20"),
                     expand = c(0,0)) +
  xlab(expression(paste("Transmission intervention efficacy % ", "(",phi,")"))) +
  ylab(expression(paste("Decolonization rate", " (", gamma[int],")"))) +
  theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank()) +
  ggtitle("Equilibrium prevalence in the community \nafter an intervention")
p.comm

pdf(file = "simul_int_eqprev_heatmap.pdf", width = 8, height = 5.5, useDingbats=F)
p.hosp
p.comm 
dev.off()


# HEATMAP OF PERCENTAGE OF CASES AVERTED  ===========================================
# both.data <- read.csv("simul_int_data.csv")

# not a perfect estimate b/c there are very slight deviations among scenarios
prop.data <- both.data %>% mutate(perc.averted = cases.averted/8223*100)

prop.data <- prop.data %>% filter(assumption != "heterogen")
prop.data$assumption <- gsub("homogen", "Homogeneous", prop.data$assumption)
prop.data$assumption <- gsub("baseline", "Risk Group", prop.data$assumption)
prop.data$assumption <- gsub("inf_adm", "Colonized Admission", prop.data$assumption)
prop.data$assumption <- gsub("inf_los", "Colonized LOS", prop.data$assumption)
prop.data$assumption <- gsub("inf_both", "Colonized LOS + Admission", prop.data$assumption)
prop.data$assumption <- as.factor(prop.data$assumption)

# smoothed with linear interpolation
p.1b <- ggplot(data = prop.data %>% filter(assumption != "Homogeneous"), 
             aes(x = beta.int*100, y = decol.rate, color = perc.averted)) +
  geom_raster(aes(fill = perc.averted), interpolate = TRUE) +
  facet_wrap(~assumption, ncol = 2, scales = "free") +
  scale_fill_viridis(option = "A", guide_legend(title = "Percentage of \nCases Averted")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100"),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     labels = c("0", "0.05", "0.1", "0.15", "0.20"),
                     expand = c(0,0)) +
  xlab(expression(paste("Transmission intervention efficacy % ", "(",phi,")"))) +
  ylab(expression(paste("Decolonization rate", " (", gamma[int],")"))) +
  theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank())
p.1b

# tiled w/o interpolation
p.2b <- ggplot(data = prop.data %>% filter(assumption != "Homogeneous"), 
              aes(x = beta.int*100, y = decol.rate, fill = perc.averted)) +
  geom_tile(aes(fill = perc.averted)) +
  facet_wrap(~assumption, ncol = 2, scales = "free") +
  scale_fill_viridis(option = "A", guide_legend(title = "Percentage of \nCases Averted")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100"),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     labels = c("0", "0.05", "0.1", "0.15", "0.20"),
                     expand = c(0,0)) +
  xlab(expression(paste("Transmission intervention efficacy % ", "(",phi,")"))) +
  ylab(expression(paste("Decolonization rate", " (", gamma[int],")"))) +
  theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank())
p.2b

# add contour lines
p.3b <- ggplot(data = prop.data %>% filter(assumption != "Homogeneous"), 
              aes(x = beta.int*100, y = decol.rate, z = perc.averted)) +
  geom_tile(aes(fill = perc.averted)) +
  stat_contour(bins = 10, color = "white") +
  facet_wrap(~assumption, ncol = 2, scales = "free") +
  scale_fill_viridis(option = "A", guide_legend(title = "Percentage of \nCases Averted")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100"),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     labels = c("0", "0.05", "0.1", "0.15", "0.20"),
                     expand = c(0,0)) +
  xlab(expression(paste("Transmission intervention efficacy % ", "(",phi,")"))) +
  ylab(expression(paste("Decolonization rate", " (", gamma[int],")"))) +
  theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank())
p.3b

# add homogeneous example
p.4b <- ggplot(data = prop.data, 
              aes(x = beta.int*100, y = decol.rate, z = perc.averted)) +
  geom_tile(aes(fill = perc.averted)) +
  stat_contour(bins = 10, color = "white") +
  facet_wrap(~assumption, ncol = 3, scales = "free") +
  scale_fill_viridis(option = "A", guide_legend(title = "Percentage of \nCases Averted")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100"),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20),
                     labels = c("0", "0.05", "0.1", "0.15", "0.20"),
                     expand = c(0,0)) +
  xlab(expression(paste("Transmission intervention efficacy % ", "(",phi,")"))) +
  ylab(expression(paste("Decolonization rate", " (", gamma[int],")"))) +
  theme_cowplot(font_size = 12) +
  theme(strip.background = element_blank()) 
p.4b


pdf(file = "simul_int_prop_cases_heatmap.pdf", width = 8, height = 5.5, useDingbats=F)
p.1b
p.2b
p.3b
p.4b
dev.off()


# HEATMAP FOR PUBLICATION  =======================================================

# save eps for submission
save_plot("heatmap.eps", p.2b,
          base_height = 4.5,
          base_width = 6.5)
