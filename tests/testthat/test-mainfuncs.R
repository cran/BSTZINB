test_that("BNB works", {

  skip_on_cran()

  data(county.adjacency)
  data(USAcities)

  # library(dplyr)
  # library(gt)
  # library(gtsummary)

  IAcities <- USAcities %>% filter(state_id=="IA")
  countyname <- unique(IAcities$county_name)
  A <- get_adj_mat(county.adjacency,countyname,c("IA"))
  m <- apply(A,1,sum)

  set.seed(091017)
  n <- nrow(A)			 # Number of spatial units
  nis <- rep(24,n) 	 # Number of individuals per county; here it's balanced -- 24 per county per year
  sid <- rep(1:n,nis)
  tid <- rep(1:nis[1],n)
  N <- length(sid)

  # library(spam)
  # library(mvtnorm)
  kappa <- 0.99			  	            # Spatial dependence parameter ~ 1 for intrinsic CAR
  covm <- matrix(c(.5,.10,.10,-.10,
                   .10,.15,.10,.10,
                   .10,.10,.5,.10,
                   -.10,.10,.10,.15),4,4)# Conditional Cov of phi1 and phi2 given phi_{-i}
  Q <- as.spam(diag(m)-kappa*A)
  covphi <- solve(Q)%x%covm			        # 3n x 3n covariance of phis
  phi <- rmvnorm(1,sigma=covphi)		    # 3n vector of spatial effects
  phitmp <- matrix(phi,ncol=4,byrow=T)  # n x 3 matrix of spatial effects

  true.phi1 <- phi1 <- phitmp[,1]-mean(phitmp[,1])   # n x 1 phi1 vector -- Centered
  true.phi2 <- phi2 <- phitmp[,2]-mean(phitmp[,2])   # n x 1 phi2 vector, etc.
  true.phi3 <- phi3 <- phitmp[,3]-mean(phitmp[,3])
  true.phi4 <- phi4 <- phitmp[,4]-mean(phitmp[,4])
  Phi1 <- rep(phi1,nis)
  Phi2 <- rep(phi2,nis)
  Phi3 <- rep(phi3,nis)
  Phi4 <- rep(phi4,nis)
  true.phi <- cbind(true.phi1,true.phi2,true.phi3,true.phi4)

  # library(boot)
  tpop <- USAcities %>% filter(state_id=="IA") %>% group_by(county_fips) %>% summarise(tpop=population) %>% .[,2] %>% unlist %>% as.numeric
  x <- rep(log(tpop),nis)
  x <- (x-mean(x))/sd(x)
  t <- tid / max(tid)
  X <- cbind(1,x)                       # Design matrix, can add additional covariates (e.g., race, age, gender)
  Xtilde <- cbind(1,x,sin(t))
  p <- ncol(Xtilde)


  # Binomial part
  true.alpha <- rep(0.1,p)
  true.alpha[1:2] <- c(-0.25,0.25)
  eta1 <- Xtilde%*%true.alpha+Phi1+Phi2*t
  pii <- inv.logit(eta1)        # 1-pr("structural zero")
  u <- rbinom(N,1,pii)          # at-risk indicator
  N1 <- sum(u)
  pstruct0 <- 1-mean(u)         # Proportion of structural zeros

  # Count Part
  true.beta <- rep(0.1,p)
  true.beta[1:2] <- c(.5,-.25)
  tmp <- u
  tmp[u==0] <- 0                       # If in structural group then set tmp to 0
  nis1 <- tapply(tmp,sid,sum)          # Number of at risk observations in each county

  eta2 <- Xtilde[u==1,]%*%true.beta+Phi3[u==1]+Phi4[u==1]*t[u==1] # Linear predictor for count part
  true.r <- 1.25                       # NB dispersion
  psi <- exp(eta2)/(1+exp(eta2))       # Prob of success
  mu <- true.r*psi/(1-psi)             # NB mean


  y <- rep(0,N)                           # Response
  y[u==1] <- rnbinom(N1,size=true.r,mu=mu)# If at risk, draw from NB

  this_dat <- data.frame(sid,tid,y,X)

  expect_snapshot({
    res0 <- BNB(y, X, A, nchain=3, niter=100, nburn=20)
  })

})

test_that("BZINB works", {

  skip_on_cran()

  data(county.adjacency)
  data(USAcities)

  # library(dplyr)
  # library(gt)
  # library(gtsummary)

  IAcities <- USAcities %>% filter(state_id=="IA")
  countyname <- unique(IAcities$county_name)
  A <- get_adj_mat(county.adjacency,countyname,c("IA"))
  m <- apply(A,1,sum)

  set.seed(091017)
  n <- nrow(A)			 # Number of spatial units
  nis <- rep(24,n) 	 # Number of individuals per county; here it's balanced -- 24 per county per year
  sid <- rep(1:n,nis)
  tid <- rep(1:nis[1],n)
  N <- length(sid)

  # library(spam)
  # library(mvtnorm)
  kappa <- 0.99			  	            # Spatial dependence parameter ~ 1 for intrinsic CAR
  covm <- matrix(c(.5,.10,.10,-.10,
                   .10,.15,.10,.10,
                   .10,.10,.5,.10,
                   -.10,.10,.10,.15),4,4)# Conditional Cov of phi1 and phi2 given phi_{-i}
  Q <- as.spam(diag(m)-kappa*A)
  covphi <- solve(Q)%x%covm			        # 3n x 3n covariance of phis
  phi <- rmvnorm(1,sigma=covphi)		    # 3n vector of spatial effects
  phitmp <- matrix(phi,ncol=4,byrow=T)  # n x 3 matrix of spatial effects

  true.phi1 <- phi1 <- phitmp[,1]-mean(phitmp[,1])   # n x 1 phi1 vector -- Centered
  true.phi2 <- phi2 <- phitmp[,2]-mean(phitmp[,2])   # n x 1 phi2 vector, etc.
  true.phi3 <- phi3 <- phitmp[,3]-mean(phitmp[,3])
  true.phi4 <- phi4 <- phitmp[,4]-mean(phitmp[,4])
  Phi1 <- rep(phi1,nis)
  Phi2 <- rep(phi2,nis)
  Phi3 <- rep(phi3,nis)
  Phi4 <- rep(phi4,nis)
  true.phi <- cbind(true.phi1,true.phi2,true.phi3,true.phi4)

  # library(boot)
  tpop <- USAcities %>% filter(state_id=="IA") %>% group_by(county_fips) %>% summarise(tpop=population) %>% .[,2] %>% unlist %>% as.numeric
  x <- rep(log(tpop),nis)
  x <- (x-mean(x))/sd(x)
  t <- tid / max(tid)
  X <- cbind(1,x)                       # Design matrix, can add additional covariates (e.g., race, age, gender)
  Xtilde <- cbind(1,x,sin(t))
  p <- ncol(Xtilde)


  # Binomial part
  true.alpha <- rep(0.1,p)
  true.alpha[1:2] <- c(-0.25,0.25)
  eta1 <- Xtilde%*%true.alpha+Phi1+Phi2*t
  pii <- inv.logit(eta1)        # 1-pr("structural zero")
  u <- rbinom(N,1,pii)          # at-risk indicator
  N1 <- sum(u)
  pstruct0 <- 1-mean(u)         # Proportion of structural zeros

  # Count Part
  true.beta <- rep(0.1,p)
  true.beta[1:2] <- c(.5,-.25)
  tmp <- u
  tmp[u==0] <- 0                       # If in structural group then set tmp to 0
  nis1 <- tapply(tmp,sid,sum)          # Number of at risk observations in each county

  eta2 <- Xtilde[u==1,]%*%true.beta+Phi3[u==1]+Phi4[u==1]*t[u==1] # Linear predictor for count part
  true.r <- 1.25                       # NB dispersion
  psi <- exp(eta2)/(1+exp(eta2))       # Prob of success
  mu <- true.r*psi/(1-psi)             # NB mean


  y <- rep(0,N)                           # Response
  y[u==1] <- rnbinom(N1,size=true.r,mu=mu)# If at risk, draw from NB

  this_dat <- data.frame(sid,tid,y,X)

  expect_snapshot({
    res1 <- BZINB(y, X, A, nchain=3, niter=100, nburn=20)
  })

})

test_that("BSTNB works", {

  skip_on_cran()

  data(county.adjacency)
  data(USAcities)

  # library(dplyr)
  # library(gt)
  # library(gtsummary)

  IAcities <- USAcities %>% filter(state_id=="IA")
  countyname <- unique(IAcities$county_name)
  A <- get_adj_mat(county.adjacency,countyname,c("IA"))
  m <- apply(A,1,sum)

  set.seed(091017)
  n <- nrow(A)			 # Number of spatial units
  nis <- rep(24,n) 	 # Number of individuals per county; here it's balanced -- 24 per county per year
  sid <- rep(1:n,nis)
  tid <- rep(1:nis[1],n)
  N <- length(sid)

  # library(spam)
  # library(mvtnorm)
  kappa <- 0.99			  	            # Spatial dependence parameter ~ 1 for intrinsic CAR
  covm <- matrix(c(.5,.10,.10,-.10,
                   .10,.15,.10,.10,
                   .10,.10,.5,.10,
                   -.10,.10,.10,.15),4,4)# Conditional Cov of phi1 and phi2 given phi_{-i}
  Q <- as.spam(diag(m)-kappa*A)
  covphi <- solve(Q)%x%covm			        # 3n x 3n covariance of phis
  phi <- rmvnorm(1,sigma=covphi)		    # 3n vector of spatial effects
  phitmp <- matrix(phi,ncol=4,byrow=T)  # n x 3 matrix of spatial effects

  true.phi1 <- phi1 <- phitmp[,1]-mean(phitmp[,1])   # n x 1 phi1 vector -- Centered
  true.phi2 <- phi2 <- phitmp[,2]-mean(phitmp[,2])   # n x 1 phi2 vector, etc.
  true.phi3 <- phi3 <- phitmp[,3]-mean(phitmp[,3])
  true.phi4 <- phi4 <- phitmp[,4]-mean(phitmp[,4])
  Phi1 <- rep(phi1,nis)
  Phi2 <- rep(phi2,nis)
  Phi3 <- rep(phi3,nis)
  Phi4 <- rep(phi4,nis)
  true.phi <- cbind(true.phi1,true.phi2,true.phi3,true.phi4)

  # library(boot)
  tpop <- USAcities %>% filter(state_id=="IA") %>% group_by(county_fips) %>% summarise(tpop=population) %>% .[,2] %>% unlist %>% as.numeric
  x <- rep(log(tpop),nis)
  x <- (x-mean(x))/sd(x)
  t <- tid / max(tid)
  X <- cbind(1,x)                       # Design matrix, can add additional covariates (e.g., race, age, gender)
  Xtilde <- cbind(1,x,sin(t))
  p <- ncol(Xtilde)


  # Binomial part
  true.alpha <- rep(0.1,p)
  true.alpha[1:2] <- c(-0.25,0.25)
  eta1 <- Xtilde%*%true.alpha+Phi1+Phi2*t
  pii <- inv.logit(eta1)        # 1-pr("structural zero")
  u <- rbinom(N,1,pii)          # at-risk indicator
  N1 <- sum(u)
  pstruct0 <- 1-mean(u)         # Proportion of structural zeros

  # Count Part
  true.beta <- rep(0.1,p)
  true.beta[1:2] <- c(.5,-.25)
  tmp <- u
  tmp[u==0] <- 0                       # If in structural group then set tmp to 0
  nis1 <- tapply(tmp,sid,sum)          # Number of at risk observations in each county

  eta2 <- Xtilde[u==1,]%*%true.beta+Phi3[u==1]+Phi4[u==1]*t[u==1] # Linear predictor for count part
  true.r <- 1.25                       # NB dispersion
  psi <- exp(eta2)/(1+exp(eta2))       # Prob of success
  mu <- true.r*psi/(1-psi)             # NB mean


  y <- rep(0,N)                           # Response
  y[u==1] <- rnbinom(N1,size=true.r,mu=mu)# If at risk, draw from NB

  this_dat <- data.frame(sid,tid,y,X)

  expect_snapshot({
    res2 <- BSTNB(y, X, A, nchain=3, niter=100, nburn=20)
  })

})

test_that("BSTZINB works", {

  skip_on_cran()

  data(county.adjacency)
  data(USAcities)

  # library(dplyr)
  # library(gt)
  # library(gtsummary)

  IAcities <- USAcities %>% filter(state_id=="IA")
  countyname <- unique(IAcities$county_name)
  A <- get_adj_mat(county.adjacency,countyname,c("IA"))
  m <- apply(A,1,sum)

  set.seed(091017)
  n <- nrow(A)			 # Number of spatial units
  nis <- rep(24,n) 	 # Number of individuals per county; here it's balanced -- 24 per county per year
  sid <- rep(1:n,nis)
  tid <- rep(1:nis[1],n)
  N <- length(sid)

  # library(spam)
  # library(mvtnorm)
  kappa <- 0.99			  	            # Spatial dependence parameter ~ 1 for intrinsic CAR
  covm <- matrix(c(.5,.10,.10,-.10,
                   .10,.15,.10,.10,
                   .10,.10,.5,.10,
                   -.10,.10,.10,.15),4,4)# Conditional Cov of phi1 and phi2 given phi_{-i}
  Q <- as.spam(diag(m)-kappa*A)
  covphi <- solve(Q)%x%covm			        # 3n x 3n covariance of phis
  phi <- rmvnorm(1,sigma=covphi)		    # 3n vector of spatial effects
  phitmp <- matrix(phi,ncol=4,byrow=T)  # n x 3 matrix of spatial effects

  true.phi1 <- phi1 <- phitmp[,1]-mean(phitmp[,1])   # n x 1 phi1 vector -- Centered
  true.phi2 <- phi2 <- phitmp[,2]-mean(phitmp[,2])   # n x 1 phi2 vector, etc.
  true.phi3 <- phi3 <- phitmp[,3]-mean(phitmp[,3])
  true.phi4 <- phi4 <- phitmp[,4]-mean(phitmp[,4])
  Phi1 <- rep(phi1,nis)
  Phi2 <- rep(phi2,nis)
  Phi3 <- rep(phi3,nis)
  Phi4 <- rep(phi4,nis)
  true.phi <- cbind(true.phi1,true.phi2,true.phi3,true.phi4)

  # library(boot)
  tpop <- USAcities %>% filter(state_id=="IA") %>% group_by(county_fips) %>% summarise(tpop=population) %>% .[,2] %>% unlist %>% as.numeric
  x <- rep(log(tpop),nis)
  x <- (x-mean(x))/sd(x)
  t <- tid / max(tid)
  X <- cbind(1,x)                       # Design matrix, can add additional covariates (e.g., race, age, gender)
  Xtilde <- cbind(1,x,sin(t))
  p <- ncol(Xtilde)


  # Binomial part
  true.alpha <- rep(0.1,p)
  true.alpha[1:2] <- c(-0.25,0.25)
  eta1 <- Xtilde%*%true.alpha+Phi1+Phi2*t
  pii <- inv.logit(eta1)        # 1-pr("structural zero")
  u <- rbinom(N,1,pii)          # at-risk indicator
  N1 <- sum(u)
  pstruct0 <- 1-mean(u)         # Proportion of structural zeros

  # Count Part
  true.beta <- rep(0.1,p)
  true.beta[1:2] <- c(.5,-.25)
  tmp <- u
  tmp[u==0] <- 0                       # If in structural group then set tmp to 0
  nis1 <- tapply(tmp,sid,sum)          # Number of at risk observations in each county

  eta2 <- Xtilde[u==1,]%*%true.beta+Phi3[u==1]+Phi4[u==1]*t[u==1] # Linear predictor for count part
  true.r <- 1.25                       # NB dispersion
  psi <- exp(eta2)/(1+exp(eta2))       # Prob of success
  mu <- true.r*psi/(1-psi)             # NB mean


  y <- rep(0,N)                           # Response
  y[u==1] <- rnbinom(N1,size=true.r,mu=mu)# If at risk, draw from NB

  this_dat <- data.frame(sid,tid,y,X)

  expect_snapshot({
    res3 <- BSTZINB(y, X, A, LinearT=TRUE, nchain=3, niter=100, nburn=20, nthin=1)
  })

})
