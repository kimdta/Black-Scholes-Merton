
#I. Merton Constraint for Call option

MertConstCall <- function(PriceOpt, S, K, r, TimeToMat){
      dimens<- length(PriceOpt)           #dimension is length of option price identified for our loop
      res   <- logical(length = dimens)   #result is a logical variable return values TRUE/FALSE
      #At each stage of the loop i= 1,2,...dimens, we check the conditions inside the loop
      for(i in c(1:dimens)){   
      # We check non-arbitrage Merton's Constraint for a Call option
      cond  <- PriceOpt[i]>= max(c(S[i]-K[i]*exp(-r[i]*TimeToMat[i]), 0))
      cond  <- cond & S[i]>= PriceOpt[i]   #union of 2 conditions is "&"
      res[i]<-cond
    }
    return(res)
}

#II. FUNCTION FOR error BS
#We build up a function to minimize error between theoretical B&S price and Market price of option which can be decomposed by 2 steps as below:

#Step 1: Black and Scholes pricing formula for European Call option:
CallBS <- function(S0, K, TimeToMat, r, Sigma){
  d1 <- (log(S0/K)+(r+1/2*Sigma^2)*TimeToMat)/(Sigma*sqrt(TimeToMat))
  d2 <- d1- Sigma*sqrt(TimeToMat)
  res<- S0*pnorm(d1)- K*exp(-r*TimeToMat)*pnorm(d2)   #using continuous compound interest regime
}

#Step 2: Compute Errors
errorBS<- function(par, S0, K, TimeToMat, r, CallPrice){
  #this function return the square of distance between quoted Call price and B&S price
  Sigma <- par  #implied volatility
  #error is a numerical object of (BSprice- Quotedprice)^2  
  error <- (CallBS(S0=S0, K=K, TimeToMat=TimeToMat, r=r, Sigma=par)- CallPrice)^2
  return(error)  
}


#III. FUNCTION FOR CALIBRATION VOLATILITY FOR BS MODEL USING MSE
#we compute Mean of errors which is an average of square distance between theoretical B&S price and market price of Call options
#Using function CallBS specified in point II

errorforCalibration<- function(par, S0, K, TimeToMat, r, CallPrice){
    #error is a numeric object of (BSprice- Quotedprice)^2  
    error <- (CallBS(S0=S0, K=K, TimeToMat=TimeToMat, r=r, Sigma=par)- CallPrice)^2
    return(mean(error))  #MSE bw B&S price and Market option price 
}

#IV. FUNCTION FOR GEOMETRIC ASIAN OPTION PRICE

AsianOptMC <- function(St0, K,  grid, Delta, Sigma, Rfree, Nsim, lev){
  #grid is a numeric object containing c(t0,t1,t2...tN)
  Ntime       <- length(grid)   #tN-t0
  GeometricSt <- rep(St0,Nsim)  #initialize object GeometricSt with i=2 we update GeometricSt= St0*St1,... i=Ntime we update to St0*St1*...*StN
  St_old      <- rep(St0,Nsim)  #initialize object St_old starting from St0, we compute trajectory of St and update St to become St_old over each iteration 
  #To generate the final realizations for the geometric prices observed on the grid, we use loop "for"
  for(i in c(2:Ntime)){ 
      #Stage 1: From St_{i-1}, we obtain St_{i} using explicit solution of GBM ; i=2 means we are computing St_{1}= St_{0}*exp(...)
      St      <- St_old*exp((Rfree- 1/2*Sigma^2)*Delta+ Sigma*sqrt(Delta)*rnorm(Nsim))
      #Stage 2: We update GeometricSt at time t_{i} by multiplying GeometricSt at time t_{i-1} with St_{i}
      GeometricSt<- GeometricSt*St  #St_{0}*St_{1}*...*St_{i} =(St_{0}*...*St_{i-1})* St_{i}
      #Stage 3: We update St to be an old price and jump into a next iteration
      St_old  <- St
  }
  #Stage 4: We determine Sample Mean
  PowerPrice <- GeometricSt^(1/Ntime)   #Ntime=N+1 in formula of final payoff
  #Discounted Final payoff under martingale measure Q
  g_n        <- exp(-Rfree*(grid[Ntime]-grid[1]))*pmax(PowerPrice-K, 0) 
  #MCPrice is an Expectation of Discounted Final payoffs is Sample Average of all values in g_n
  MCPrice    <- mean(g_n)
  #Stage 5: We determine MC bounds for 95% level
  alpha      <- (1-lev)/2      #2.5%  
  # Determine quantiles for the bounds
  q_alpha_lb <- qnorm(alpha)   # -1.96     
  q_alpha_ub <- -q_alpha_lb    # +1.96
  # Sample variance for the Discounted Final payoff
  Samplevar  <- var(g_n) 
  # Upper/ lower bounds of MCPrice
  LB <- MCPrice + q_alpha_lb*1/sqrt(Nsim)*sqrt(Samplevar)
  UB <- MCPrice + q_alpha_ub*1/sqrt(Nsim)*sqrt(Samplevar)
  
  return(list(MCPrice= MCPrice, LowerB= LB, UpperB= UB, DeltaB= abs(UB-LB)))
}




