###############################################################################
# Discrete to continous
# Approach: Rolling window
# Period data
###############################################################################

rm(list = ls())
library(tidyverse)
library(data.table)
library(zoo)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("Mx.RData") # load data for females


Roller <- function(Rx, Year, nrol, fn = sum, skip){
  rol <-rollapply(Rx, nrol, fn, by = skip)
  Year.rol <- c(seq(Year[1],Year[1] + length(rol) -1, by = 1 ))
  out<- data.frame(Year=Year.rol, rol)
  return(out)}


Mx <- data.table(Mx)

rDx <- as.data.frame(Mx[, Roller(Rx= Dx, Year = Year, nrol = 10, skip = 1),
                        by = list(Age, country, Sex)])  
rNx <- as.data.frame(Mx[, Roller(Rx= Nx, Year = Year, nrol = 10, skip = 1),
                        by = list(Age, country, Sex)])

names(rDx)[5] <- "Dx"
names(rNx)[5] <- "Nx"

rMx <- merge(rDx, rNx, by = c("Age", "country", "Sex", "Year"))

rMx$Mx <- rMx$Dx / rMx$Nx
rMx <- data.table(arrange(rMx, country, Year, Age))


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

# Comparison with rho calculated with raw data
rho <- rMx[, get.rho(mx = Mx, Year= Year, Age= Age ), by = list(country, Sex)]
rho.mx <- Mx[, get.rho(mx = Mx, Year= Year, Age= Age ), by = list(country, Sex)]

ggplot()+
  geom_line(data = subset(rho.mx, Age == 10 & Sex == "Females"
                          & country == "Japan"), aes(Year, rho))+
    geom_line(data = subset(rho, Age == 10
                            & Sex == "Females"
                            & country == "Japan"), aes(Year, rho), colour = "red", size = 1.3)+
  theme_minimal()+
  coord_cartesian(xlim = c(1950,2020))+
  theme(panel.grid.minor = element_blank())+
  ylab("Rates of mortality improvement at age 80")


