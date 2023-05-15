source("Functions.R")

#Upload dataset
library(readxl)
DatasetGroups1_4 <- read_excel("DatasetGroups1-4.xlsx")
attach(DatasetGroups1_4)

#Store variables into new objects
S0   <- DatasetGroups1_4$Underlying
K    <- DatasetGroups1_4$Strike
Call <- DatasetGroups1_4$CallPrice
TimeToMat <- DatasetGroups1_4$TimeToMat/ 360  #by convention a year= 360days
r   <- DatasetGroups1_4$rfree
length(Call)    #20 Call options

###QUESTION I. CLEAN DATASET

#Step 1.Check the no-arbitrage bounds (Merton constraint)
CheckedCall<- MertConstCall(PriceOpt= Call, S= S0, K= K, 
                  r= r, TimeToMat= TimeToMat)
CheckedCall # All Quoted Options satisfy Merton Constraint.

#Step 2. Check Monotonicity (by graphical approach)
plot(K, Call)   #behavior of Call price is non-increasing with respect to K

#Step 3. (In case CheckCall return some value of FALSE) we remove option which does not satisfy Merton constraint (remove options returned value of FALSE)
PriceOptCLeaned <- Call[CheckedCall]
KCleaned <- K[CheckedCall]
CleanedDataset <-DatasetGroups1_4[CheckedCall,]   #extract all TRUE elements

#Step 4. Compute percentage of options in cleaned data set:
length(KCleaned)/length(K)   #100% cleaned so we can use original dataset


###QUESTION II. REPRESENT VOLATILITY SMILE
#Volatility Smile: For each quoted price, we obtain only 1 implicit volatility for all maturities

par0 <- 0.1   #initial guess for value of implied volatility
ImpliedVol <- numeric(length=length(Call))  #initialize numeric object contains implied volatilities
ValueofObj <- numeric(length=length(Call))  #object contains values of objective function

#For each quoted option, we perform a minimization problem with respect to sigma using an objective function of errorBS
for (i in c(1:length(Call))){
  res <- optim(par = par0, fn= errorBS, method = "L-BFGS-B",
               lower= 0.001, S0=S0[i], K=K[i], TimeToMat= TimeToMat[i],
               r= r[i], CallPrice= Call[i])
  if(res$convergence==0){ #if convergence=0, we have successful completion   
    cat("\n", "successful convergence")
  }else{   #otherwise checking the position of option and its value of error
    cat("\n Check the prob numb.", c(i,res$value))
  }
  ImpliedVol[i] <- res$par     #solution of optimization problem
  ValueofObj[i] <- res$value   #values of Objective function
}

ImpliedVol  #implied volatility ranging from 0.413 to 0.369
ValueofObj  #values of Objective function are all close to 0

plot(x= K, y= ImpliedVol, main= "Volatility Smile for Call Option")
#Implied volatility is monotonically decreasing with respect to K

###QUESTION III. CALIBRATE THE VOLATILITY 

#1. Calibrate the volatility for BS model using Mean Square Error

#Calibration: we get an unique sigma from quoted options which is optimal using MSE, passing function errorforCalibration
res<- optim(par=0.10, fn=errorforCalibration, lower=0.01,
           method ="L-BFGS-B",
           S0=S0, K=K, TimeToMat=TimeToMat, 
           r=r, CallPrice=Call)

CalibratedImpVol <- res$par; CalibratedImpVol  
res$value         #error using BS is quite high
res$convergence   #convergence=0

#2. Pricing Butterfly strategy using calibrated volatility

#Define strike prices
K1 <- 4165        #1 Long position in a Call with Strike price K1
K2 <- 4175        #1 Long position in a Call with Strike price K2
Kbar <- 4170      #2 Short position in a Call with Strike price Kbar (also recognize that Kbar=(K1+K2)/2)

#Pricing each components in Butterfly strategy using B&S with same underlying price and calibrated implied volatility
S0  <- S0[1]          # take first value as a represented value of Underlying price
r   <- r[1]           # take first value as a represented value of risk free rate
TimeToMat <- 25/360   #transform into yearly basis

CallBS1 <- CallBS(S0= S0, K= K1, TimeToMat= TimeToMat, r=r, Sigma= CalibratedImpVol); CallBS1
CallBS2 <- CallBS(S0= S0, K= K2, TimeToMat= TimeToMat, r=r, Sigma= CalibratedImpVol); CallBS2
CallBS3 <- CallBS(S0= S0, K= Kbar, TimeToMat= TimeToMat, r=r, Sigma= CalibratedImpVol); CallBS3

#Non-arbitrage price at time t of Butterfly strategy can be constructed as:
Butterflyprice<- 1*CallBS1 + 1*CallBS2 - 2*CallBS3
Butterflyprice

#3. Probability of having strictly positive Final payoff from Butterfly strategy

#Step 1. Final payoff of Butterfly strategy by graphical representation:

FinalPayoff<-function(S_T, K1, K2, Kbar){
  finalpayoff <- pmax(S_T-K1,0) + pmax(S_T-K2, 0) - 2*pmax(S_T- Kbar, 0)
  plot(S_T, finalpayoff, type="l")
}
Butterfly<- FinalPayoff(S_T =seq(4150,4190, by=1), K1=K1, K2= K2, Kbar= Kbar)

#We observe that Strictly positive final payoff of Butterfly is obtained when: K1 < S_T < K2

#Step 2. Applying conditional probability of a stochastic process of Underlying price over interval (K1,K2) given info at time t0
#S_T follows a GBM with drift of risk free rate; 

Sigma <- CalibratedImpVol   #using Calibrated implied volatility
Timetomat<- 25/360          #transfrom Timetomat into yearly basis
St0   <- S0[1]              # take first value as a represented value of Underlying price
rfree <- r[1]               # take first value as a represented value of risk free rate

Prob <- function(St0, K1, K2, Sigma, Timetomat){
  d2 <- (log(K2/St0)-(rfree-0.5*Sigma^2)*Timetomat)/(Sigma*sqrt(Timetomat))
  d1 <- (log(K1/St0)-(rfree-0.5*Sigma^2)*Timetomat)/(Sigma*sqrt(Timetomat))
  Pr <- pnorm(d2) - pnorm(d1)
  return(Pr)
}
#Probability of having strictly positive payoff using Calibrated implied volatility is 0.9267%
Prob(St0= St0, K1= K1, K2= K2, Sigma= CalibratedImpVol, Timetomat = Timetomat)

#Step 3. Use Vectorize() to see the behavior of this probability as Sigma changes

ProbVect <- Vectorize(Prob, vectorize.args="Sigma")
#Define a sequence of sigma to see behavior in probability of strictly positive payoff
sigma   <- seq(0.01, 0.9, by= 0.001)  
probs   <- ProbVect(St0= St0, K1= K1, K2= K2, Sigma= sigma, Timetomat = Timetomat)

#Plot a graph of probability of having strictly positive final payoff in volatility
plot(sigma, probs, type="l")
abline(v= CalibratedImpVol, col= "red") #threshold of Calibrated implied volatility   

###IV. COMPUTE MONTE CARLO PRICE FOR GEOMETRIC ASIAN OPTION

#We simulate 1000 random samples to compute Monte Carlo price for Geometric Asian Option 
Nsim      <- 1000       #number of simulations
TimeToMat <- 30/360     #transform TimeToMat into yearly basis
Delta     <- 5/360      #Each sub-interval will be equal to: t_{i} - t_{i-1}
grid      <- seq(0, TimeToMat, by= Delta) #equally spaced time grid: t_{1}, t_{2},...,t_{i-1}, t_{i},...,t_{N}
GeomK     <- 4170       #strike price of Geometric Asian option

#Geometric Asian Option price using function AsianOptMC with 95% confidence level
#We use set.seed() to fix the generated random sample, each time we run we will obtain the same result
set.seed(123)  
AsianOptMC(St0= St0, K=GeomK,  grid= grid, Delta = Delta,
           Sigma= CalibratedImpVol, Rfree= rfree, Nsim= Nsim, lev= 0.95)

#We use Vectorize() to see the behavior of Asian MC Price as the Nsim increases
set.seed(123)
AsianOptVect <- Vectorize(AsianOptMC, vectorize.args="Nsim")
AsianOptVect(St0=St0, K=GeomK,  grid= grid, Delta= Delta,
             Sigma=CalibratedImpVol, Rfree=rfree, Nsim= 10^(3:8), lev = 0.95)

# As number of simulation (Nsim) increases from 1000 to 10^8, MCPrice converge to 102.547 
#and length of interval between LB and UP gets smaller

