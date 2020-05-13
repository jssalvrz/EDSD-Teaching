###############################################################################
# Estimation of Gamma Gompertz and Gamma Gompertz via the mode
# Period data
###############################################################################
rm(list = ls())
library(tidyverse)
library(data.table)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("Dat.RData")
source("GG_FUN.R") # Load Gamma-Gompertz Functions via the mode

# Fit Gamma-Gompertz for one population -----------------------------------

Japan <- subset(Dat, country == "Japan" & Year == 2010 & Sex == "Females" & Age>=80)

fit.Japan <- GG(D = Japan$Dx,E = Japan$Nx) # Fit Gamma Gompertz 
fit.Japan$Age <- fit.Japan$Age + 79.5 # Rescale ages

# Parameters
fit.Japan[1,c(5,6,7,8,9)]

# Theoretical value of the plateau b/gamma
fit.Japan$b[1] / fit.Japan$gamma[1] 

# Plot the estimated hazard vs the raw death rates
ggplot()+
  geom_point(data = Japan, aes(Age,Mx), size = 1, colour = "grey20")+
  geom_line(data = fit.Japan, aes(Age, hx), colour = "red")+
  geom_ribbon(data = fit.Japan, aes(Age,hx, ymin = hx.low, ymax = hx.upper),
              alpha = 0.2, fill = "red")+
  scale_y_continuous(trans = "log", breaks = c(0.01, 0.05,0,1, 0.2,0.5,0.7), expand = c(0,0))+
  scale_x_continuous(breaks = seq(70,105, by = 5), expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(xlim = c(80,105), ylim = c(min(Japan$Mx), 0.7))+
  theme(panel.grid.minor = element_blank())+
  ylab("Mortality hazard (log)")

# Fit Gamma-Gompertz for all countries and years --------------------------

Age.o <- 80 # starting age to fit GG. What is the "best" age to start the analysis?
Dat <- data.table(subset(Dat, Age>= Age.o & Age <=105)) # Select ages

fit <- data.frame(Dat[,optim(lnL,x=0:(105-Age.o),
                             D=Dx, # Death counts
                             E=Nx, # Exposures
                             control=list(fnscale=-1), # Maximization
                             ageInterval=1, # Length of the interval
                             par=c(0.5,0.1,0.3))$par, # Initial parameters
                             by=list(Year,Sex,country)])

# Tidy up parameters
fit$par <- c("a", "b", "gamma") # Rename parameters
fit <- fit %>% spread(par, V1)
fit <- arrange(fit, Sex, country, Year)

# Plot parameters of the Gamma-Gompertz
# a vs b
ggplot(fit)+
  geom_point(aes(a,b, colour = Sex))+
  #facet_wrap(~country)+
  theme_minimal()+
  ggtitle("Parameters Gamma-Gompertz (a,b)")

# b vs gamma
ggplot(fit)+
  geom_point(aes(gamma,b, colour = Sex))+
  #facet_wrap(~country)+
  theme_minimal()+
  ggtitle("Parameters Gamma-Gompertz (a,gamma)")

# Rate of ageing over time
ggplot(fit)+
  geom_line(aes(Year, b, colour = Sex))+
  facet_wrap(~country)+
  scale_y_continuous(breaks = seq(0.08,0.20, by =0.02), expand = c(0,0))+
  scale_x_continuous(breaks = seq(1950,2020, by = 10), expand = c(0,0))+
  coord_cartesian(ylim = c(0.08,0.2))+
  theme_minimal()+
  ggtitle("Rate of ageing over time, Gamma-Gompertz")

# Fit Gamma-Gompertz via the Mode -----------------------------------------

Mfit <- data.frame(Dat[,optim(MlnL,x=0:(105-Age.o),
                             D=Dx,
                             E=Nx,
                             control=list(fnscale=-1),
                             ageInterval=1,
                             par=c(10,0.1,0.3))$par, by=list(Year,Sex, country)])

# Tidy up the parameters
Mfit$par <- c("M", "b", "gamma") # Rename parameters
Mfit <- Mfit %>% spread(par, V1)
Mfit <- arrange(Mfit, Sex,country, Year)
Mfit$M <- Mfit$M + Age.o


# Compare both models -----------------------------------------------------


# Plot parameters of the Gamma-Gompertz
# b vs gamma
ggplot(Mfit)+
  geom_point(aes(b,gamma, colour = Sex))+
  #facet_wrap(~country)+
  theme_minimal()+
  coord_cartesian(ylim = c(0,0.3))+
  ggtitle("Parameters Gamma-Gompertz via the Mode (b, gamma)")

# Rate of ageing over time
ggplot(Mfit)+
  geom_line(aes(Year, b, colour = Sex))+
  facet_wrap(~country)+
  scale_y_continuous(breaks = seq(0.08,0.20, by =0.02), expand = c(0,0))+
  scale_x_continuous(breaks = seq(1950,2020, by = 10), expand = c(0,0))+
  coord_cartesian(ylim = c(0.08,0.2))+
  theme_minimal()+
  ggtitle("Rate of ageing over time, Gamma-Gompertz")

# Comparison of rate of ageing
ggplot()+
  geom_line(data = subset(fit, Sex == "Females"),aes(Year, b), colour = "forestgreen")+
  geom_line(data = subset(Mfit, Sex == "Females"),aes(Year, b))+
  facet_wrap(~country)+
  scale_y_continuous(breaks = seq(0.08,0.20, by =0.02), expand = c(0,0))+
  scale_x_continuous(breaks = seq(1950,2020, by = 10), expand = c(0,0))+
  coord_cartesian(ylim = c(0.08,0.2))+
  theme_minimal()+
  ggtitle("Rate of ageing over time, GG vs GGM. Females")

# Modal age at death
ggplot(subset(Mfit, Sex == "Females"))+
  geom_line(aes(Year, M, colour = country))+
 # facet_wrap(~country)+
  scale_y_continuous(breaks = seq(70,100, by =2), expand = c(0,0))+
  scale_x_continuous(breaks = seq(1950,2020, by = 10), expand = c(0,0))+
  coord_cartesian(ylim = c(78,92))+
  theme_minimal()+
  ylab("Modal age at death (years)")+
  ggtitle("Modal age at death, Females")
  
