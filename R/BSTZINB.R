#' @name BSTZINB
#' @title Fit a Bayesian Spatiotemporal Zero Inflated Negative Binomial model
#'
#' @description
#' Generate posterior samples for the parameters in a Bayesian Spatiotemporal Zero Inflated Negative Binomial Model
#'
#' @usage BSTZINB(y,X,A,LinearT = TRUE,
#'             nchain=3,niter=100,nburn=20,nthin=1)
#'
#' @param y vector of counts, must be non-negative
#' @param X matrix of covariates, numeric
#' @param A adjacency matrix, numeric
#' @param LinearT logical, whether to fit a linear or non-linear temporal trend
#' @param nchain positive integer, number of MCMC chains to be run
#' @param niter positive integer, number of iterations in each chain
#' @param nburn non-negative integer, number of iterations to be discarded as burn-in samples
#' @param nthin positive integer, thinning interval
#'
#' @importFrom stats cov
#' @importFrom stats dnbinom
#' @importFrom stats rbinom
#' @importFrom stats reorder
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats spline
#' @importFrom stats var
#' @import dplyr
#' @import BayesLogit
#' @import spam
#' @import MCMCpack
#' @import msm
#' @import splines
#' @import boot
#' @importFrom matrixcalc is.positive.definite
#'
#' @return list of posterior samples of the parameters of the model
#'
#' @examples
#' data(simdat)
#' y <- simdat$y
#' X <- cbind(simdat$V1,simdat$x)
#' data(county.adjacency)
#' data(USAcities)
#' IAcities <- subset(USAcities,state_id=="IA")
#' countyname <- unique(IAcities$county_name)
#' A <- get_adj_mat(county.adjacency,countyname,c("IA"))
#' \donttest{
#' res3 <- BSTZINB(y, X, A, LinearT=TRUE, nchain=2, niter=100, nburn=20, nthin=1)
#' }
#'
#' @export
BSTZINB <- function(y, X, A, LinearT = TRUE, nchain=3, niter=100, nburn=20, nthin=1){

  ## Run the necessary checks
  if(!is.vector(y)){stop("y must be a vector")}
  if(!is.matrix(X)){stop("X must be a matrix")}
  if(!is.matrix(A)){stop("A must be a matrix")}
  if(nchain < 1){stop("nchain must be a positive integer")}
  if(niter < 1){stop("niter must be a positive integer")}
  if(nburn < 0){stop("nburn must be a non-negative integer")}
  if(nthin < 1){stop("nthin must be a positive integer")}
  y <- as.numeric(y)
  if(min(y,na.rm=T)<0){stop("y must be non-negative")}
  if(!is.numeric(X)){stop("X must be numeric")}
  if(!is.numeric(A)){stop("A must be numeric")}
  if(!is.logical(LinearT)){stop("LinearT must be logical")}

  N <- length(y)
  n <- nrow(A)			    # Number of spatial units
  nt <- N/n
  nis <- rep(nt,n) 		# Number of individuals per county; here it's balanced
  # Note: may need to lower proposal variance, s, below as n_i increases
  sid <- rep(1:n,nis)
  tid <- rep(1:nis[1],n)
  t <- tid / max(tid)
  if(is.null(colnames(X))){colnames(X) <- c("intercept", paste0("X",1:(ncol(X)-1)))}
  if(LinearT==T){
    Xtilde <- cbind(X,t)
  }else{
    Tbs <- bs(t,df=7)
    colnames(Tbs) <- paste0("t",colnames(Tbs))
    Xtilde <- cbind(X,Tbs)
  }
  p <- ncol(Xtilde)

  ##########
  # Priors #
  ##########
  alpha0 <- beta0 <- rep(0,p)
  T0a <- diag(.01,p)
  T0b <- diag(.01,p)         # Uniform or Gamma(0.01,0.01) prior for r depending on MH or Gibbs
  s <- 0.0003                # Proposal variance  -- NOTE: may need to lower this as n_i increases
  kappa <- 0.999999
  Q <- as.spam(diag(apply(A,1,sum))) - kappa*as.spam(A)

  ############
  # Num Sims #
  ############
  lastit <- (niter - nburn)/nthin	# Last stored value

  ############
  # Store #
  ############
  Beta <- Alpha <- array(0,c(lastit,p,nchain))
  colnames(Beta) <- colnames(Alpha) <- colnames(Xtilde)
  R <- matrix(0,lastit,nchain)
  Sigphi <- array(0,c(lastit,16,nchain))
  PHI1 <- PHI2 <- PHI3 <- PHI4 <- array(0,c(lastit,n,nchain))
  I <- Eta1 <- Eta2 <- array(0,c(lastit,N,nchain))

  for(chain in 1:nchain){

    #########
    # Inits #
    #########

    beta <- alpha <- rnorm(p)
    phi_init <- rmvnorm(1,sigma=diag(.1,16*n))	  # Random effects
    phi_init <- matrix(phi_init,ncol=16,byrow=T)  # n x 3 matrix of spatial effects
    phi1 <- phi_init[,1]
    phi2 <- phi_init[,2]
    phi3 <- phi_init[,3]
    phi4 <- phi_init[,4]
    phi1 <- phi1-mean(phi1)
    phi2 <- phi2-mean(phi2)
    phi3 <- phi3-mean(phi3)
    phi4 <- phi4-mean(phi4)
    Phi1 <- rep(phi1,nis)
    Phi2 <- rep(phi2,nis)
    Phi3 <- rep(phi3,nis)
    Phi4 <- rep(phi4,nis)

    phimat <- cbind(phi1,phi2,phi3,phi4)
    Sigmaphi <- cov(phimat)

    r <- 1
    Acc <- 0
    y1 <- rep(0,N)             # At risk indicator (this is W in paper)
    y1[y > 0] <- 1               # If y>0, then at risk w.p. 1
    N0 <- length(y[y == 0])      # Number of observed 0's
    q <- rep(.5,N)             # 1-p=1/(1+exp(X%*%alpha)), used for updating y1

    ########
    # MCMC #
    ########

    for (i in 1:niter){

      # Update alpha
      mu <- Xtilde%*%alpha+Phi1+Phi2*t
      w <- rpg(N,1,mu)
      z <- (y1-0.5)/w
      v <- solve(crossprod(sqrt(w)*Xtilde)+T0a)
      m <- v%*%(T0a%*%alpha0+t(sqrt(w)*Xtilde)%*%(sqrt(w)*(z-Phi1-Phi2*t)))
      alpha <- c(rmvnorm(1,m,v))

      # Update phi1
      priorprec <- as.numeric(1/(Sigmaphi[1,1]-Sigmaphi[1,-1]%*%solve(Sigmaphi[-1,-1])%*%Sigmaphi[-1,1]))*Q # Prior Prec of phi1|phi2,phi3,phi4
      priormean <- diag(n)%x%(Sigmaphi[1,-1]%*%solve(Sigmaphi[-1,-1]))%*%c(t(phimat[,-1]))      # Prior mean of phi1|phi2,phi3,phi4
      prec <- priorprec+as.spam(diag(tapply(w,sid,sum),n,n))
      m <- c(priorprec%*%priormean)+tapply(w*(z-Xtilde%*%alpha-Phi2*t),sid,sum)
      if(is.positive.definite(prec%>%as.matrix)) phi1 <- rmvnorm.canonical(1, m, prec)[1,]

      # Center
      phi1 <- phi1-mean(phi1)
      Phi1 <- rep(phi1,nis)

      # Update phi2
      priorprec <- as.numeric(1/(Sigmaphi[2,2]-Sigmaphi[2,-2]%*%solve(Sigmaphi[-2,-2])%*%Sigmaphi[-2,2]))*Q # Prior Prec of phi2|phi1,phi3,phi4
      priormean <- diag(n)%x%(Sigmaphi[2,-2]%*%solve(Sigmaphi[-2,-2]))%*%c(t(phimat[,-2]))      # Prior mean of phi2|phi1,phi3,phi4
      prec <- priorprec+as.spam(diag(tapply(w*t^2,sid,sum),n,n))
      m <- c(priorprec%*%priormean)+tapply(t*w*(z-Xtilde%*%alpha-Phi1),sid,sum)
      if(is.positive.definite(prec%>%as.matrix)) phi2 <- rmvnorm.canonical(1, m, prec)[1,]

      # Center
      phi2 <- phi2-mean(phi2)
      Phi2 <- rep(phi2,nis)

      # Update r
      rnew <- rtnorm(1,r,sqrt(s),lower=0)       # Treat r as continuous
      ratio <- sum(dnbinom(y[y1==1],rnew,q[y1==1],log=T)) - sum(dnbinom(y[y1==1],r,q[y1==1],log=T))+
        dtnorm(r,rnew,sqrt(s),0,log=T) - dtnorm(rnew,r,sqrt(s),0,log=T)   # Uniform Prior for R
      # Proposal not symmetric
      if (log(runif(1)) < ratio) {
        r <- rnew
        Acc <- Acc+1
      }

      # Update at-risk indicator y1 (W in paper)
      eta1 <- Xtilde%*%alpha + Phi1 + Phi2*t
      eta2 <- Xtilde%*%beta + Phi3 + Phi4*t              # Use all n observations
      pii <- pmax(0.01,pmin(0.99,inv.logit(eta1)))  # at-risk probability
      q <- pmax(0.01,pmin(0.99,1/(1+exp(eta2))))                      # Pr(y=0|y1=1)
      theta <- pii*(q^r)/(pii*(q^r)+1-pii)         # Conditional prob that y1=1 given y=0 -- i.e. Pr(chance zero|observed zero)
      y1[y==0] <- rbinom(N0,1,theta[y==0])      # If y=0, then draw a "chance zero" w.p. theta, otherwise y1=1
      N1 <- sum(y1)
      nis1 <- tapply(y1,sid,sum)

      # Update beta
      eta <- Xtilde[y1==1,]%*%beta + Phi3[y1==1] + Phi4[y1==1]*t[y1==1]
      w <- rpg(N1,y[y1==1]+r,eta)                               # Polya weights
      z <- (y[y1==1]-r)/(2*w)
      v <- solve(crossprod(Xtilde[y1==1,]*sqrt(w))+T0b)
      m <- v%*%(T0b%*%beta0+t(sqrt(w)*Xtilde[y1==1,])%*%(sqrt(w)*(z-Phi3[y1==1]-Phi4[y1==1]*t[y1==1])))
      beta <- c(rmvnorm(1,m,v))

      # Update phi3
      n1 <- length(nis1[nis1>0])

      priorprec <- as.numeric(1/(Sigmaphi[3,3]-Sigmaphi[3,-3]%*%solve(Sigmaphi[-3,-3])%*%Sigmaphi[-3,3]))*Q # Prior Prec of phi3|phi1,phi2,phi4
      priormean <- diag(n)%x%(Sigmaphi[3,-3]%*%solve(Sigmaphi[-3,-3]))%*%c(t(phimat[,-3]))      # Prior mean of phi3|phi1,phi2,phi4
      prec <- priorprec+as.spam(diag(tapply(w,sid[y1==1],sum),n,n))
      tmp <- rep(0,n)                       # Account empty blocks
      tmp[nis1>0] <- tapply(w*(z-Xtilde[y1==1,]%*%beta-Phi4[y1==1]*t[y1==1]),sid[y1==1],sum)
      m <- c(priorprec%*%priormean)+tmp
      if(is.positive.definite(prec%>%as.matrix)) phi3 <- rmvnorm.canonical(1, m, prec)[1,]

      # Center
      phi3 <- phi3-mean(phi3)
      Phi3 <- rep(phi3,nis)

      # Update phi4
      priorprec <- as.numeric(1/(Sigmaphi[4,4]-Sigmaphi[4,-4]%*%solve(Sigmaphi[-4,-4])%*%Sigmaphi[-4,4]))*Q # Prior Prec of phi4|phi1,phi2,phi3
      priormean <- diag(n)%x%(Sigmaphi[4,-4]%*%solve(Sigmaphi[-4,-4]))%*%c(t(phimat[,-4]))      # Prior mean of phi4|phi1,phi2,phi3

      prec <- priorprec+as.spam(diag(tapply(w*t[y1==1]^2,sid[y1==1],sum),n,n))
      tmp <- rep(0,n)                       # Account for empty counties
      tmp[nis1>0] <- tapply(t[y1==1]*w*(z-Xtilde[y1==1,]%*%beta-Phi3[y1==1]),sid[y1==1],sum)
      m <- c(priorprec%*%priormean)+tmp
      if(is.positive.definite(prec%>%as.matrix)) phi4 <- rmvnorm.canonical(1, m, prec)[1,]

      # Center
      phi4 <- phi4-mean(phi4)
      Phi4 <- rep(phi4,nis)

      # Update r2 using Gibbs as in Dadaneh et al and Zhou and Carin
      # Update latent counts, l
      # l <- rep(0,N1)
      # ytmp <- y[y1==1]
      # for(j in 1:N1) l[j] <- sum(rbinom(ytmp[j],1,round(r/(r+1:ytmp[j]-1),6))) # Could try to avoid loop; rounding avoids numerical stability

      # Update r from conjugate gamma distribution given l and psi
      # psi <- exp(eta2[y1==1])/(1+exp(eta2[y1==1]))
      # r2 <- rgamma(1,0.01+sum(l),0.01-sum(log(1-psi)))


      # Update Sigma.phi
      phimat <- cbind(phi1,phi2,phi3,phi4)
      try({
        Sigmaphi <- riwish(4+n-1,diag(4)+t(phimat)%*%Q%*%phimat)
      })

      # Store
      if (i> nburn & i%%nthin==0) {
        j <- (i-nburn)/nthin
        Alpha[j,,chain] <- alpha
        Beta[j,,chain] <- beta
        R[j,chain] <- r
        # R2[j,chain] <- r2
        Sigphi[j,,chain] <- c(Sigmaphi)
        PHI1[j,,chain] <- phi1
        PHI2[j,,chain] <- phi2
        PHI3[j,,chain] <- phi3
        PHI4[j,,chain] <- phi4
        Eta1[j,,chain] <- eta1
        Eta2[j,,chain] <- eta2
        I[j,,chain] <- y1
      }

      # if (i%%10==0) print(paste(chain, "/", nchain,"chain | ",round(i/niter*100,2),"% completed"))
      # if (i%%10==0) print(paste(chain, "/", nchain,"chain | ",round(i/niter*100,2),"% completed |","Test:",conv.test(R[,chain])))

    }
  }

  list_params = list(Alpha=Alpha, Beta=Beta, R=R, Sigphi=Sigphi,
                     PHI1=PHI1, PHI2=PHI2, PHI3=PHI3, PHI4=PHI4,
                     Eta1=Eta1, Eta2=Eta2, I=I)
  class(list_params) <- "DCMB"
  return(list_params)

}
