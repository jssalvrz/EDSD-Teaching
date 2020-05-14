###############################################################################
# Discrete to continous
# Approach: Smooth data in two dimensions (age and time)
# Period data
###############################################################################

rm(list = ls())
library(tidyverse)
library(data.table)
library(splines)
library(MortalitySmooth)
library(pracma)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("Mx.RData") # load data for females


# Convert dataframes into matrices ----------------------------------------

# We use Japanese females as an example
Dx <- subset(Mx[,-c(6,7)],country == "Japan" & Sex== "Females")
Nx <- subset(Mx[,-c(5,7)],country == "Japan" & Sex == "Females")


# Convert data frames into matrices
Dx <- as.data.frame(spread(Dx[,-c(1,2)],Year, Dx))
rownames(Dx) <- Dx$Age
Dx[Dx == 0] <- 0.00001 # To avoid NAs
Dx <- as.matrix(Dx[,-1])

Nx <- as.data.frame(spread(Nx[,-c(1,2)],Year,Nx))
rownames(Nx) <- Nx$Age
Nx[Nx == 0] <- 0.00001 # To avoid NAs
Nx <- as.matrix(Nx[,-1])

# Smooth death counts and calculate death rates ---------------------------
# We need to define the ages and years that will be part of the smoothing
Ages  <- as.integer(rownames(Dx)) 
Years <- as.integer(colnames(Dx))

# Smooth deaths
fit <- Mort2Dsmooth(x = Ages, y = Years, 
                      Z = Dx, offset = log(Nx))

sDx <- fit$fitted.values # These are the smoothed death counts
smx <- as.data.frame(sDx/Nx) # Smoothed death rates
smx$Age <- as.integer(rownames(smx))

# Convert matrix to dataframe again
smx <- gather(smx, Year, Mx,-Age)
smx$Year <- as.integer(smx$Year)

# Compare results

ggplot(subset(Mx, Sex == "Females" & Year == 2010 & country == "Japan"))+
  geom_line(aes(Age, Mx))+
  geom_line(data = subset(smx, Year == 2010), aes(Age, Mx), colour = "red")+
  scale_y_continuous(trans = "log", breaks = c(0.01, 0.05,0,1, 0.2,0.5,0.7), expand = c(0,0))+
  scale_x_continuous(breaks = seq(0,110, by = 10), expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(xlim = c(0,110), ylim = c(min(Mx$Mx), 0.7))+
  theme(panel.grid.minor = element_blank())+
  ylab("Mortality hazard (log)")


ggplot(subset(Mx, Sex == "Females" & Age == 80 & country == "Japan"))+
  geom_line(aes(Year, Mx))+
  geom_line(data = subset(smx, Age == 80), aes(Year, Mx), colour = "red")+
  scale_y_continuous(trans = "log", breaks = c(0.01, 0.05,0,1, 0.2,0.5,0.7), expand = c(0,0))+
  scale_x_continuous(breaks = seq(1950,2020, by = 10), expand = c(0,0))+
  theme_minimal()+
  coord_cartesian(xlim = c(1950,2020))+
  theme(panel.grid.minor = element_blank())+
  ylab("Mortality hazard (log)")

# Continous death rates ---------------------------------------------------

# Function to interpolate 
interp <- function(fx, Age, s.age = 0){
  ispl <- interpSpline( fx ~ Age, bSpline = TRUE )
  x    <- predict(ispl, seq( s.age, 110, 0.01 ))$x
  y    <- predict(ispl, seq( s.age, 110, 0.01 ))$y
  out  <- data.frame(Age= x, hx = y)
  return(out)}

smx <- data.table(smx)

# Calculation for all years
smx.inter<- smx[,interp(fx =  Mx, Age = Age, s.age = 0),  by = list(Year)]


# Rates of mortality improvement ------------------------------------------

# This function is to calculate rathes of mortality improvement (rho)
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

rho <- smx.inter[, get.rho(mx = hx, Year= Year, Age= Age )]

# Comparison with rho calculated with unsmoothed data
Mxj <- data.table(subset(Mx, country == "Japan" & Sex == "Females" ))
rho.mx <- Mxj[, get.rho(mx = Mx, Year= Year, Age= Age )]

ggplot()+
  geom_line(data = subset(rho.mx, Age == 80), aes(Year, rho))+
  geom_line(data = subset(rho, Age == 80), aes(Year, rho), colour = "red")+
  theme_minimal()+
  coord_cartesian(xlim = c(1950,2020))+
  theme(panel.grid.minor = element_blank())+
  ylab("Rates of mortality improvement at age 80")


# Modal age at death and percentiles --------------------------------------

# Function to calculate survival, density and cumulative hazard functions
calc.Sx <- function(Age, hx){
  
  Hx <- cumtrapz(Age,hx)
  Sx <- exp(-Hx)
  Fx <- 1-Sx
  fx <- hx * Sx
  dat <- data.frame(Age,hx,Hx,Sx,Fx,fx)
  
  return(dat)}


Sx <- smx.inter[, calc.Sx(Age=Age,hx=hx), by = list(Year)]

# Function to calculate the modal age at death
get.mode <- function(Age, fx){
  
  r<- data.frame(Age, fx)
  M <- r$Age[r$fx == max(r$fx)]
  return(M)
  
}

Mode <- Sx[, get.mode(Age= Age, fx = fx), by =list(Year)]
names(Mode) <- c("Year", "Mode")

# Modal age at death over time
ggplot(Mode)+
  geom_line(aes(Year, Mode))+
  scale_y_continuous(breaks = seq(0,110, by = 5), expand = c(0,0))+
  scale_x_continuous(breaks = seq(1950,2020, by = 10), expand = c(0,0))+
  #coord_cartesian(xlim = c(1950,2020), ylim = c(75,95) )+
  theme_minimal()
