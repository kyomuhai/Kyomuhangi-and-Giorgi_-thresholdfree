#******************************************************************************************** 
#******************************************************************************************** 
####BACKGROUND 
# 
# This script contains syntax used for the analysis of malaria serology data as described in the paper: 
# "A threshold-free approach with age-dependency for estimating malaria seroprevalence", as specified for the M2 approach. 

# The antibody measurements used in this analysis are PfAMA OD values obtained from ELISA, 
# however the methods are applicable to any malaria antigen type, and continuous antibody measurement. 

# Throughout the script we indicate where the code may need to change depending on the dataset/antibody type under analysis. 
# Explanations of the functions and operations used in the sytax are provided, 
# and further details of statistical/mathematical principles of this analysis can be found in the paper. 

# To request access to the dataset used, please contact Gillian Stresman (Gillian.Stresman@lshtm.ac.uk) or Chris Drakeley (Chris.Drakeley@lshtm.ac.uk) at LSHTM.  

#******************************************************************************************** 
#*#******************************************************************************************** 



rm(list=ls())
setwd("~/ama analysis") #set the working directory accordingly

#load the required dataset 
load("amadata.RData") 
# this dataset should contain the continuous antibody measurements and age of individuals
# ensure that each observation is not missing any of these two variables. 


#********************************************************************************************      
#Generate required vectors and objects
#********************************************************************************************  
y <- log(data$ama_norm)
a <- data$age
ind.age <- list()
max.age <- max(a)

#********************************************************************************************      
#Create a function to implement equations 3, 8 and 10. 
#******************************************************************************************** 

# For this analysis, we selected a linear spline with a knot at age 10 for equation 8, based on results from figure 3(a) (see paper for details). 
# This may need to be modified depending on the dataset and/or antigen being analysed. 
# We also select age as a logit linear predictor in equation 10, to account for age dependency. 
# Other functional forms for these predictors can be explored (see paper for details)

#********************************************* 
##Linear predictors as decribed in equations 8 and 10
#*********************************************  
#equation 8: 
formula.mean.y <- ~ age + I((age>10)*(age-10))
#equation 10: 
formula.mix.prob <- ~ age



#*****************************************************************************************   
## Combining equations 3, 8 and 10 in the Empirical Model (as illustrated in Fig 1) 
#*****************************************************************************************    

lr.estim.mixture.emp <- function(formula.mean.y,
                                 formula.mix.prob,messages=TRUE,
                                 start.theta=NULL,
                                 return.hessian=FALSE)  {
  a.max <- max(a)
  #matrix for mu(a) in equation 8:
  D.t.y <- model.matrix(formula.mean.y,data=data.frame(age=a))
  #matrix for p(a) in equation 10:
  D.t.mix.prob <- model.matrix(formula.mix.prob,data=data.frame(age=a))
  
  time.covar <- 0:(a.max-1)
  #number of params for equation 8: 
  p.mean.y <- ncol(D.t.y)
  #number of params for equation 10: 
  p.mix.prob <- ncol(D.t.mix.prob)
  #initialize starting parameters at 0 for all parameters:
  if(is.null(start.theta)) start.theta <- c(rep(0,p.mean.y+p.mix.prob+3)) 
  #accounting for truncation of data 
  threshold <- tapply(y,a,max)
  threshold.data <- sapply(a,function(i)
    threshold[i])
  #log likelihood function for the empirical model (illustrated in Figure 1):
  log.lik <- function(theta) {
    #equation 8 params: 
    beta <- theta[1:p.mean.y]
    #equation 10 params: 
    beta.tilde <- theta[(p.mean.y+1):(p.mean.y+p.mix.prob)]
    #equation 3 params:
    delta <- 1+exp(theta[p.mean.y+p.mix.prob+1])
    sigma2.mix1 <- exp(theta[p.mean.y+p.mix.prob+2])
    sigma2.mix2 <- exp(theta[p.mean.y+p.mix.prob+3])
    #equation 10:
    pr.a.mix <- 1/(1+exp(-D.t.mix.prob%*%beta.tilde))
    #equation 3 params:
    mu.a.mix1 <- exp(D.t.y%*%beta)
    mu.a.mix2 <- delta*mu.a.mix1
    #natural log transformation of the mean and sd in equation 3:
    mean.ln.mix1 <- log(mu.a.mix1/sqrt(1+sigma2.mix1/(mu.a.mix1^2)))
    sd.ln.mix1 <- sqrt(log(1+sigma2.mix1/(mu.a.mix1^2)))
    mean.ln.mix2 <- log(mu.a.mix2/sqrt(1+sigma2.mix2/(mu.a.mix2^2)))
    sd.ln.mix2 <- sqrt(log(1+sigma2.mix2/(mu.a.mix2^2)))
    #accounting for truncation in the distribution function:
    trunc.prob <- pr.a.mix*pnorm(threshold.data,
                                 mean=mean.ln.mix2,
                                 sd=sd.ln.mix2,log=FALSE)+
      (1-pr.a.mix)*pnorm(threshold.data,
                         mean=mean.ln.mix1,
                         sd=sd.ln.mix1,log=FALSE)
    #probability density:
    num <- (pr.a.mix*dnorm(y,
                           mean=mean.ln.mix2,
                           sd=sd.ln.mix2,log=FALSE)+
              (1-pr.a.mix)*dnorm(y,
                                 mean=mean.ln.mix1,
                                 sd=sd.ln.mix1,log=FALSE))
    #the likelihood, accounting for truncation: 
    llik <- sum(log(num)-log(trunc.prob))
  }
  #optimization of the log likelihood: 
  estim <- nlminb(start.theta,
                  function(x) -log.lik(x),
                  control=list(trace=1*messages))  
  
  library(numDeriv)
  if(return.hessian) estim$H <- hessian(log.lik,estim$par)
  #extracting the mixture model parameters (see table 2)
  estim$regression.mean.antibody <- estim$par[1:p.mean.y] 
  estim$regression.mixing.prob <- estim$par[(p.mean.y+1):(p.mean.y+p.mix.prob)] 
  names(estim$regression.mean.antibody) <- colnames(D.t.y) 
  names(estim$regression.mixing.prob) <- colnames(D.t.mix.prob)
  estim$mean.fact <-  1+exp(estim$par[p.mean.y+p.mix.prob+1])
  estim$sigma2.mix1 <-  exp(estim$par[p.mean.y+p.mix.prob+2]) 
  estim$sigma2.mix2 <-  exp(estim$par[p.mean.y+p.mix.prob+3]) 
  estim$formula.mean.y <- formula.mean.y 
  estim$formula.mix.prob <- formula.mix.prob 
  estim$D.t.y <- D.t.y 
  estim$D.t.mix.prob <- D.t.mix.prob 
  #calculate the aic
  estim$aic <- 2*length(estim$par)+2*estim$objective 
  return(estim)
}



#********************************************************************************************      
#Fit the Empirical Model 
#********************************************************************************************  

estim.emp <- lr.estim.mixture.emp(formula.mean.y,
                                  formula.mix.prob,
                                  return.hessian = TRUE)



#********************************************************************************************      
#Extract CIs from the estim.emp object 
#********************************************************************************************  
# use the Hessian matrix to obtain CIs for the parameters. 
# If you are maximising a likelihood, then the covariance matrix of the estimates 
# is (asymptotically) the inverse of the negative of the Hessian. 
# The standard errors are the square roots of the diagonal elements of the covariance
fisher_info_emp <- solve(-estim.emp$H) 
prop_sigma_emp<- sqrt(diag(fisher_info_emp))

emp.params  <- data.frame (estim.emp$par)
emp.params.prop <- data.frame(prop_sigma_emp)
names(emp.params)[1] <- "par"
names(emp.params.prop)[1] <- "prop_sigma_emp"

emp.interval <- merge(emp.params, emp.params.prop, by=0, all=TRUE)
emp.interval$lower_emp <- emp.interval$par-1.96*emp.interval$prop_sigma_emp
emp.interval$upper_emp <- emp.interval$par+1.96*emp.interval$prop_sigma_emp
emp.interval <-  emp.interval %>% 
  select(-c(prop_sigma_emp, "Row.names"))
emp.interval <- t(emp.interval)
emp.interval <- as.data.frame(emp.interval)
# Follow the proceedure from the lr.estim.mixture.emp function when extracting the parameter estimates. 
#Note the transformation of delta, sigma2.mix1 and sigma2.mix2 from the log scale. These need to be transformed in the same way below 
emp.interval$V6 <- 1+exp(emp.interval$V6)
emp.interval$V7 <- exp(emp.interval$V7)
emp.interval$V8 <- exp(emp.interval$V8)
emp.interval
# double check that these values match the parameters (the values and order) specified in the estim.emp object. 



#********************************************************************************************      
# Plot the mixture distribution 
#********************************************************************************************  
age.dist <- function(estim.emp,age,plot.legend=FALSE) {
  D.t.y <- model.matrix(estim.emp$formula.mean.y,data=data.frame(age=age))
  
  D.t.mix.prob <- model.matrix(estim.emp$formula.mix.prob,data=data.frame(age=age))
  
  
  mu.age.mix1 <- exp(D.t.y%*%estim.emp$regression.mean.antibody)
  mu.age.mix2 <- mu.age.mix1*estim.emp$mean.fact
  
  mean.g.mix1.emp = log(mu.age.mix1/sqrt(1+estim.emp$sigma2.mix1/(mu.age.mix1^2)))
  sd.g.mix1.emp = sqrt(log(1+estim.emp$sigma2.mix1/(mu.age.mix1^2)))
  
  mean.g.mix2.emp =log(mu.age.mix2/sqrt(1+estim.emp$sigma2.mix2/(mu.age.mix2^2)))
  sd.g.mix2.emp =sqrt(log(1+estim.emp$sigma2.mix2/(mu.age.mix2^2)))
  
  threshold <- tapply(y,a,max)
  threshold.data <- threshold[age]
  
  pr.t.mix.emp <- as.numeric(1/(1+exp(-D.t.mix.prob%*%estim.emp$regression.mixing.prob)))
  cut.off <- (mean.g.mix1.emp+(3*sd.g.mix1.emp))
  
  par(mar=c(5.1,5,4.1,2.5))
  hist(y[a==age], xlab="log OD", probability = TRUE, ylab="",
       main=paste("Age ",age),breaks=10, cex.lab=2, cex.axis=1.5, cex.main=2, 
       cex.sub=1.5, cex.axis=2,xlim=c(-8,3), ylim=c(0,0.5))
  
  
  f.emp <-  function(x) exp(log(pr.t.mix.emp*dnorm(x,mean.g.mix2.emp,sd.g.mix2.emp)+
                                  (1-pr.t.mix.emp)*dnorm(x,mean.g.mix1.emp,sd.g.mix1.emp))-
                              log(pr.t.mix.emp*pnorm(threshold.data,mean.g.mix2.emp,sd.g.mix2.emp)+
                                    (1-pr.t.mix.emp)*
                                    pnorm(threshold.data,mean.g.mix1.emp,sd.g.mix1.emp)))
  f.emp <- Vectorize(f.emp,"x")
  curve(f.emp,col=4,add = TRUE, lwd=3,lty="solid",xlim = c(-10,threshold.data))
  if(plot.legend) legend(-1,0.5,c("Unified","Empirical"),col=c(2,4),lty=c("solid","dashed"),lwd=2.5,cex=1.5)
  abline(v=cut.off, col="red",lty="dashed",lwd =2)
  text(x = cut.off+1.5, y = 0.45, label=round(cut.off, digits = 2),col="red",cex=2.5)
}

age.dist(estim.emp,age=1,plot.legend = F)
age.dist(estim.emp,age=2,plot.legend = F)
age.dist(estim.emp,age=3,plot.legend = F)
age.dist(estim.emp,age=4,plot.legend = F)
age.dist(estim.emp,age=5,plot.legend = F)
age.dist(estim.emp,age=6,plot.legend = F)
age.dist(estim.emp,age=7,plot.legend = F)
age.dist(estim.emp,age=8,plot.legend = F)
age.dist(estim.emp,age=9,plot.legend = F)
age.dist(estim.emp,age=10,plot.legend = F)
age.dist(estim.emp,age=11,plot.legend = F)
age.dist(estim.emp,age=12,plot.legend = F)
age.dist(estim.emp,age=13,plot.legend = F)
age.dist(estim.emp,age=14,plot.legend = F)
age.dist(estim.emp,age=15,plot.legend = F)
age.dist(estim.emp,age=16,plot.legend = F)



#********************************************************************************************      
# Extract posterior probabilities  using estim.emp
#********************************************************************************************  
# Create a function which extracts posterior probabilities  using estim.emp. Much of this code is copied from the lr.estim.mixture.emp function     
post.prob <- function(estim, y.dens,age) {
  D.age <- model.matrix(estim$formula.mix.prob,data=data.frame(age=age))
  pr.t.mix <- 1/(1+exp(-as.numeric(D.age%*%estim$regression.mixing.prob)))
  
  D.t.y <- model.matrix(estim$formula.mean.y,data=data.frame(age=age))
  mu.t.mix1 <- exp(D.t.y%*%estim$regression.mean.antibody) 
  mu.t.mix2 <- estim$mean.fact*mu.t.mix1 
  
  mean.ln.mix1 <- log(mu.t.mix1/sqrt(1+estim$sigma2.mix1/(mu.t.mix1^2))) 
  sd.ln.mix1 <- sqrt(log(1+estim$sigma2.mix1/(mu.t.mix1^2))) 
  mean.ln.mix2 <- log(mu.t.mix2/sqrt(1+estim$sigma2.mix2/(mu.t.mix2^2))) 
  sd.ln.mix2 <- sqrt(log(1+estim$sigma2.mix2/(mu.t.mix2^2))) 
  
  post.prob <- (pr.t.mix*dnorm(y.dens,
                               mean=mean.ln.mix2,
                               sd=sd.ln.mix2,log=FALSE)/
                  ((pr.t.mix*dnorm(y.dens,
                                   mean=mean.ln.mix2,
                                   sd=sd.ln.mix2,log=FALSE)+
                      (1-pr.t.mix)*dnorm(y.dens,
                                         mean=mean.ln.mix1,
                                         sd=sd.ln.mix1,log=FALSE))))
  return(post.prob)
}

post.prob <- Vectorize(post.prob,c("y.dens","age"))
post.prob.obs <- post.prob(estim=estim.emp,y.dens = y,age=a)

##Vizualize posterior probabilities by age
plot.post.prob.by.age <- function(age) {
  hist(post.prob.obs[a==age],main=paste("Age",age), ylim = c(0,400), xlim = c(0,1),
       breaks=50, cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=1.5, cex.axis=2, xlab = "Posterior probabilities")
}
plot.post.prob.by.age(1)
plot.post.prob.by.age(2)
plot.post.prob.by.age(3)
plot.post.prob.by.age(4)
plot.post.prob.by.age(5)
plot.post.prob.by.age(6)
plot.post.prob.by.age(7)
plot.post.prob.by.age(8)
plot.post.prob.by.age(9)
plot.post.prob.by.age(10)
plot.post.prob.by.age(11)
plot.post.prob.by.age(12)
plot.post.prob.by.age(13)
plot.post.prob.by.age(14)
plot.post.prob.by.age(15)
plot.post.prob.by.age(16)



#************************************************************************************************************      
# Estimate seroprevalence using Monte Carlo methods (i,e, binomial sampling with posterior probabilities)  
#************************************************************************************************************  
# Create dataset containing the log transformed OD values (y) , age, age and posterior probabilities
post.df <- data.frame(y,a,post.prob.obs) 
range(a)

# Organise this dataset by age
library(tidyverse)
post.df.split <- post.df %>%
  group_split(a) 

# Create function for implementing binomial sampling which will generate 10000 realizations for each observation
bin.sampling <- function(x, nsims = 10000) {
  samples <- matrix(NA, nrow = nrow(x), ncol = nsims)
  for (i in 1:nrow(x)) {
    samples[i, ] <- rbinom(n = nsims, size = 1, prob = x$post.prob.obs)
    #print(i)
  }
  return(samples)
}
samples.dist.list <- lapply(post.df.split, bin.sampling)

# Create a dataset which combines post.df.split with the 10000 realizations
binom.list <- list()
for (i in 1:length(post.df.split)) {
  binom.list[[i]] <- cbind(post.df.split[[i]], samples.dist.list[[i]])
}

##generate the seroprevalence distribution for each age 
sero.prev.dist <- lapply(binom.list, function(x) colMeans(x[4:ncol(binom.list[[1]])]))

##Plot seroprevalence
sero.prev.plot.list <- list(mu= rep(NA,length(sero.prev.dist)),
                            q25 = rep(NA,length(sero.prev.dist)),
                            q75 = rep(NA,length(sero.prev.dist)))
for(i in 1:length(sero.prev.dist)){
  sero.prev.plot.list$mu[i] <- mean(sero.prev.dist[[i]])
  sero.prev.plot.list$q25[i] <- quantile(sero.prev.dist[[i]], prob=c(0.25))
  sero.prev.plot.list$q75[i] <- quantile(sero.prev.dist[[i]], prob=c(0.75))
}
plot(sero.prev.plot.list$mu, ylab="P", xlab="age", col="blue",pch=16, ylim=c(0,1),
     cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=1.5, cex.axis=2)
require(plotrix)
plotCI(sero.prev.plot.list$mu,li=sero.prev.plot.list$q25,ui=sero.prev.plot.list$q75, 
       lwd=1, slty=2,ylab="P", xlab="age", col="dark blue",pch=16, ylim=c(0,1),
       cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5,xaxt='n')
axis(side=1, at=seq(0,16, 2),cex.axis=1.5)



#********************************************************************************************     
# Application: Input the seroprevalence estimates in RCM 
#********************************************************************************************  


#***************************************************************************************
#Create a function to implement equation 7. 
#***************************************************************************************
# Set the regression for mu(a) to be the same as in the estim.emp object, 
# (i.e. the mixture model) In this case the formula is equation 8
formula.mean.y <- estim.emp$formula.mean.y

# Set the regression for the seroconversion rate (lambda). 
# In this analysis, we use an intercept only model for lambda, i.e there are no covariates modifying lambda.
# However, in future analysis, covariates such as altitude can be added to this formula if the interest is to measure and account for these covariates. 
# For example, in malaria epidemiology altitude can be an important predictor for transmission intensity. 
# Whether or not to include additional covariates can be tested and results compared using an appropriate index such as AIC. 
formula.log.lambda <- ~ 1

rcm.model <- function(y.bin,
                      formula.mean.y,
                      formula.log.lambda,
                      omega.fixed,
                      messages=TRUE,
                      start.theta=NULL,
                      return.hessian=FALSE) {
  a.max <- max(a)
  #matrix for lambda:
  D.lambda <- model.matrix(formula.log.lambda,data=data)
  #number or params for formula.log.lambda equation:
  p.lambda <- ncol(D.lambda)
  #initialize starting parameters at 0 for formula.log.lambda equation:
  if(is.null(start.theta)) start.theta <- c(rep(0,p.lambda)) 
  #omega is fixed
  omega <- omega.fixed
  #log likelihood function for the RCM
  log.lik <- function(theta) {
    gamma <- theta[1:p.lambda]
    lambda <- exp(D.lambda%*%gamma)
    # equation 7:
    pr.t.mix <- exp((log(lambda)-log((lambda+omega))))*(1-exp(-(lambda+omega)*a))
    llik <- sum(dbinom(y.bin,size=1,prob=pr.t.mix,log=TRUE)) 
    llik 
  }
  #optimization of the log likelihood:
  estim <-  nlminb(start.theta,
                   function(x) -log.lik(x),
                   control=list(trace=1*messages))  
  
  library(numDeriv)
  if(return.hessian) estim$H <- hessian(log.lik,estim$par)
  #extracting the RCM parameters (see table 2 - in this case it's just lambda)
  estim$regression.log.lambda <- estim$par[1:p.lambda]
  names(estim$regression.log.lambda) <- colnames(D.lambda)
  estim$omega <- omega.fixed
  estim$formula.log.lambda <- formula.log.lambda
  estim$D.lambda <- D.lambda
  #calculate the aic
  estim$aic <- 2*length(estim$par)+2*estim$objective
  return(estim)
}

#Pass the rcm function through the 10000 realizations of the binomial sampling
binom.samples.rcm <- rbindlist(binom.list)
a <- binom.samples.rcm$a
samples <- subset(binom.samples.rcm, select=-c(y,a,post.prob.obs))
lambda.val <- list()
log.lik.hat <- list()
aic.hat <- list()
col_names <- colnames(samples)

for(i in colnames(samples)) {
  y.bin <- samples[[i]]
  estim.rcm <- rcm.model(y.bin=y.bin,
                         formula.mean.y=formula.mean.y,
                         formula.log.lambda=formula.log.lambda,
                         omega.fixed = 0.01)
  lambda.hat <- exp(estim.rcm$regression.log.lambda)
  log.lik <- estim.rcm$objective
  aic <- estim.rcm$aic
  
  lambda.val[i] <- lambda.hat
  log.lik.hat[i] <- -log.lik
  aic.hat[i] <- aic
}

#extract the distribution of lambda estimates (n=10,000)
lambda.dist<- unlist(lambda.val)


#plot the lambda distribution
library(latex2exp)
par(mar=c(5.1, 4.9, 4, 2.1), mgp=c(3, 0.5, 0),las=1)
hist(lambda.dist, pch = 19, xlab=expression(lambda), main="distribution of lambda",
     cex.lab=2,cex.axis=1.5, cex=1.5, breaks=40, ylim=c(0,2000),xaxt='n' )
axis(side=1, at=seq(-3,2,0.1),cex.axis=1.5)
abline(v=mean(lambda.dist), col="blue",lty="dashed",lwd =3)
abline(v=quantile(lambda.dist,0.025), col="red",lty="dashed",lwd =2)
abline(v=quantile(lambda.dist,0.975), col="red",lty="dashed",lwd =2)


##Generate 95% CIs
lambda <- mean(lambda.dist)
lower <- quantile(lambda.dist,0.025)
upper <- quantile(lambda.dist,0.975)
interval <- data.frame(value=lambda, lower=lower, upper=upper)
interval

#Plot seroprevalence curve from RCM
seroprev.b <- as.data.frame(sero.prev.plot.list)
library(data.table)
setDT(seroprev.b, keep.rownames = TRUE)[]

compute.pr.a <- function(interval,seroprev.b) {
  lambda_v <- interval$value
  lambda_l <- interval$lower
  lambda_u <- interval$upper
  omega <- estim.rcm$omega
  pr.mu <- seroprev.b$mu
  a <- unique(seroprev.b$rn)
  
  plot(pr.mu, ylab="", xlab="age", col="blue",pch=16, ylim=c(0,1),
       cex=1, cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=1.5, cex.axis=2)
  
  f <- function(a) exp((log(lambda_v)-log((lambda_v+omega))))*(1-exp(-(lambda_v+omega)*a))
  f <- Vectorize(f,"a")
  
  f_l <- function(a) exp((log(lambda_l)-log((lambda_l+omega))))*(1-exp(-(lambda_l+omega)*a))
  f_l <- Vectorize(f_l,"a")
  
  f_u <- function(a) exp((log(lambda_u)-log((lambda_u+omega))))*(1-exp(-(lambda_u+omega)*a))
  f_u <- Vectorize(f_u,"a")
  
  
  curve(f,col="#AB1779",add = TRUE, lwd=2)
  curve(f_l,col="#AB1779",add = TRUE, lwd=1, lty=3)
  curve(f_u,col="#AB1779",add = TRUE, lwd=1, lty=3)
}

compute.pr.a(interval,seroprev.b)


save.image("ama threshold free.RData")

