# Simulation of a sample path for a Brownian Motion
set.seed(123) 
N <- 100 # number of end-points of the grid including T 
Time <- 1 # length of the interval [0,T] in time units 
Delta <- Time/N # time increment 
W <- numeric(N+1) # initialization of the vector W 
t <- seq(0, Time, length=N+1) 
for(i in 2:(N+1)){ 
  W[i] <- W[i-1] + rnorm(1) * sqrt(Delta) 
} 
plot(t, W, type = "l", main = "Brownian Motion", 
     ylim = c(-1, 1), cex.main = 0.8)

# Simulation of M sample paths of a BM with scale and drift
# Inputs:
t0=0
FinalT=1
M=100
XIntial = rep(0,M)
N=1000
sigma=0.1
mu =0.01
#Code
set.seed(123)
DeltaT <- (FinalT-t0)/N
BMwithDrift<- matrix(NA,M,N+1)
BMwithDrift[,1]<-XIntial
for(t in c(2:(N+1))){
  BMwithDrift[,t]<-BMwithDrift[,t-1]+mu*(DeltaT)+sigma*sqrt(DeltaT)*rnorm(M)
}
matplot(x = seq(t0,FinalT,by=DeltaT),
        y=t(BMwithDrift),
        type="l")
source("FunctionsSC.R")
set.seed(123)
x11()
BMsp<-SamPathBM(X0=0, t0=0, FinalT=1, M=100,N=1000,
                    sigma=sigma, mu = mu)

set.seed(123)
x11()
GBMsp<-SamPathGBM(X0=1, t0=0, FinalT=1, M=1000,N=100,
                sigma=sigma, mu = mu)
