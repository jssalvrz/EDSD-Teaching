###############################################################################
# Discrete to continous
# Approach: Simulation of lifespans using exponential distribution
#           with piece-wise constant hazard
# Period data
###############################################################################

rm(list = ls())
library(tidyverse)
library(data.table)
library(msm)
library(pracma)
library(splines)
library(mltools)
library(Rmisc)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("Mx.RData") # load data
Mx <- data.table(Mx)

# Simulate lifespans ------------------------------------------------------
# Function to simulate a great number of lifespans
sim <- function(nsim, Mx, s.age){
  l <- length(Mx)-1
  lt <- as.data.frame(rpexp(n=nsim,rate=Mx,0:l))# ideally 100000
  names(lt) <- "lt"
  lt$lt <- lt$lt + s.age
  return(lt)}

# Simulate lifespans, it takes a lot of time so I recommend to run it and safe them
sims <- Mx[, sim(nsim = 100000, Mx = Mx, s.age = 0), by =list(country, Sex, Year)]

# Calculate probability functions -----------------------------------------

# Calculate hx, Sx, Fx and fx

calc.Sx <- function(lt, Qx){ # Qx is the starting age of calculations
  dat     <- data.frame(lt, Qx)
  dat$id  <- ifelse(dat$lt>= dat$Qx, 1,0)
  dat <- subset(dat, id == 1)
  dat.cdf <- empirical_cdf(dat$lt, 
                           ubounds=seq(round(min(dat$lt), 2),
                                       round(max(dat$lt), 2)-0.01,.01)) # Calculate CDF
  
  Age <- dat.cdf$UpperBound
  Fx  <- dat.cdf$CDF
  Sx  <- 1- Fx
  
  dens <- density(dat$lt)
  fx.raw <- dens$y
  Age.raw <- dens$x
  s.age <- min(Age)
  l.age <- max(Age)
  
  ispl <- interpSpline( fx.raw ~ Age.raw, bSpline = TRUE )
  x<- predict(ispl, seq( s.age, l.age, 0.01 ))$x
  y<- predict(ispl, seq( s.age, l.age, 0.01 ))$y
  fx <- pmax(0,y)
  hx <- fx/Sx
  
  out <- data.frame(Age, Fx, Sx, fx, hx)
  
  return(out)}

Sx <- sims[,calc.Sx(lt = lt, Qx=79.98), # I'm starting at that age but it can be at any age
          by = list(country, Sex, Year)]


# Rates of mortality improvement ------------------------------------------

get.rho <- function(mx, Year, Age){ 
  MX<- data.frame(mx, Year, Age)
  MX <- MX %>% spread(Year, mx)
  mx.2   <- MX[ ,-c(1,2)]
  mx.1   <- MX[ ,-c(1,ncol(MX))]
  rho <- cbind(Age=MX$Age,-log(mx.2/mx.1)) # rho populaton by age
  rho <- rho %>% gather(Year, rho, -Age )
  rho$Year <- as.integer(rho$Year)
  return(data.frame(rho))
} 

rho <- Sx[, get.rho(mx = hx, Year= Year, Age= Age ), by = list(country, Sex)]
rho.mx <- Mx[, get.rho(mx = Mx, Year= Year, Age= Age ), by = list(country, Sex)]

# Comparison of results
ggplot()+
  geom_line(data = subset(rho.mx, Age == 80 &
                            Sex == "Females" & country == "Japan"), aes(Year, rho))+
   geom_line(data = subset(rho, Age == 80 & 
                             Sex == "Females" & country == "Japan"), aes(Year, rho), colour = "red")+
  theme_minimal()+
  coord_cartesian(xlim = c(1950,2020))+
  theme(panel.grid.minor = element_blank())+
  ylab("Rates of mortality improvement at age 80")



# Calculate life expectancy, e-dagger, entropy ----------------------------

JPN <- subset(sims, country == "Japan" & Sex == "Females" & Year == 2010)

CI(JPN$lt, 0.95) # Life expectancy and confidence intervals

# Now let's calculate it for all the countries and years

# Calculate e-dagger, H and ex
calc.H <- function(Age, Sx){
  ex <- trapz(Age, Sx) # life expectancy; trapz is used to integrate
  e.dagger <- -trapz(Age, Sx * log(Sx))
  H <- e.dagger / ex
  out <- data.frame(ex, e.dagger, H)
  return(out)
}

# Summary measures at age 79.98. They can be calculated at any age
H <- Sx[, calc.H(Age= Age, Sx = Sx), by = list(country, Sex, Year)]


ggplot(H)+
  geom_line(aes(Year, ex, colour = country))+
  facet_wrap(~Sex)+
  theme_minimal()+
  coord_cartesian(xlim = c(1950,2020))+
  theme(panel.grid.minor = element_blank())+
  ylab("Life expectancy at age 79.89")
  
ggplot(H)+
  geom_line(aes(Year, H, colour = country))+
  facet_wrap(~Sex)+
  theme_minimal()+
  coord_cartesian(xlim = c(1950,2020))+
  theme(panel.grid.minor = element_blank())+
  ylab("Entropy at age 79.89")

ggplot(H)+
  geom_point(aes(e.dagger,ex, colour = country))+
  facet_wrap(~Sex)+
  theme_minimal()+
  theme(panel.grid.minor = element_blank())+
  xlab("e-dagger at age 79.89")+
  ylab("Life expectancy at age 79.89")


ggplot(H)+
  geom_point(aes(H,ex, colour = country))+
  facet_wrap(~Sex)+
  theme_minimal()+
  theme(panel.grid.minor = element_blank())+
  xlab("Entropy at age 79.89")+
  ylab("Life expectancy at age 79.89")

