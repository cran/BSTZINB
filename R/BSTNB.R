#' @name BSTNB
#' @title Fit a Bayesian Spatiotemporal Negative Binomial model
#'
#' @description
#' Generate posterior samples for the parameters in a Bayesian Spatiotemporal Negative Binomial Model
#'
#' @usage BSTNB(y,X,A,
#'             nchain=3,niter=100,nburn=20,nthin=1)
#'
#' @param y vector of counts, must be non-negative
#' @param X matrix of covariates, numeric
#' @param A adjacency matrix, numeric
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
#' res2 <- BSTNB(y, X, A, nchain=2, niter=100, nburn=20, nthin=1)
#' }
#'
#' @export
BSTNB <- function(y, X, A, nchain=3, niter=100, nburn=20, nthin=1){

  ## Run the necessary checks
  if(!is.vector(y)){stop("y must be a vector")}
  if(!is.matrix(X)){stop("X must be a matrix")}
  if(!is.matrix(A)){stop("A must be a matrix")}
  if(nchain < 1){stop("nchain must be a positive integer")}
  if(niter < 1){stop("nsim must be a positive integer")}
  if(nburn < 0){stop("nburn must be a non-negative integer")}
  if(nthin < 1){stop("nthin must be a positive integer")}
  y <- as.numeric(y)
  if(min(y,na.rm=T)<0){stop("y must be non-negative")}
  if(!is.numeric(X)){stop("X must be numeric")}
  if(!is.numeric(A)){stop("A must be numeric")}

  N <- length(y)
  n <- nrow(A)			    # Number of spatial units
  nt <- N/n
  nis <- rep(nt,n) 		# Number of individuals per county; here it's balanced -- 50 per county per year
  # Note: may need to lower proposal variance, s, below as n_i increases
  sid <- rep(1:n,nis)
  tid <- rep(1:nis[1],n)
  t <- tid / max(tid)
  N <- length(sid) 		  # Total number of observations
  p <- ncol(X)
  if(is.null(colnames(X))){colnames(X) <- c("intercept", paste0("X",1:(ncol(X)-1)))}

  ##########
  # Priors #
  ##########
  beta0 <- rep(0,p)
  T0a <- diag(.01,p)
  T0b <- diag(.01,p)         # Uniform or Gamma(0.01,0.01) prior for r depending on MH or Gibbs
  s <- 0.0003                # Proposal variance  -- NOTE: may need to lower this as n_i increases
  kappa <- 0.999999
  Q <- as.spam(diag(apply(A,1,sum)))-kappa*as.spam(A)

  ############
  # Num Sims #
  ############
  lastit <- (niter-nburn)/nthin	# Last stored value

  ############
  # Store #
  ############
  Beta <- array(0,c(lastit,p,nchain))
  colnames(Beta) <- colnames(X)
  R <- matrix(0,lastit,nchain)
  Sigphi <- array(0,c(lastit,4,nchain))
  PHI3 <- PHI4 <- array(0,c(lastit,n,nchain))
  Eta <- array(0,c(lastit,N,nchain))

  for(chain in 1:nchain){

    #########
    # Inits #
    #########
    beta<-rnorm(p)
    phi_init <- spam::rmvnorm(1,sigma=diag(.1,2*n))	  # Random effects
    phi_init <- matrix(phi_init,ncol=2,byrow=T)  # n x 3 matrix of spatial effects
    phi3 <- phi_init[,1]
    phi4 <- phi_init[,2]
    phi3 <- phi3-mean(phi3)
    phi4 <- phi4-mean(phi4)
    Phi3 <- rep(phi3,nis)
    Phi4 <- rep(phi4,nis)

    phimat <- cbind(phi3,phi4)
    Sigmaphi <- cov(phimat)
    r <- 1
    Acc <- 0
    N0 <- length(y[y==0])      # Number of observed 0's
    q <- rep(.5,N)             # 1-p=1/(1+exp(X%*%alpha)), used for updating y1


    ########
    # MCMC #
    ########

    for (i in 1:niter){

      # Update r
      rnew <- rtnorm(1,r,sqrt(s),lower=0)       # Treat r as continuous
      ratio <- sum(dnbinom(y,rnew,q,log=T))-sum(dnbinom(y,r,q,log=T))+
        dtnorm(r,rnew,sqrt(s),0,log=T) - dtnorm(rnew,r,sqrt(s),0,log=T)   # Uniform Prior for R
      # Proposal not symmetric
      if (log(runif(1)) < ratio) {
        r <- rnew
        Acc <- Acc+1
      }

      # Update r2 using Gibbs as in Dadaneh et al and Zhou and Carin
      # Update latent counts, l
      # l <- rep(0,N)
      # ytmp <- y
      # for(j in 1:N) l[j] <- sum(rbinom(ytmp[j],1,round(r/(r+1:ytmp[j]-1),6))) # Could try to avoid loop; rounding avoids numerical stability

      # Update r from conjugate gamma distribution given l and psi
      # psi<-exp(eta2)/(1+exp(eta2))
      # r2<-rgamma(1,0.01+sum(l),0.01-sum(log(1-psi)))

      # Update beta
      eta <- X%*%beta+Phi3+Phi4*t
      w <- rpg(N,y+r,eta)                               # Polya weights
      z <- (y-r)/(2*w)
      v <- solve(crossprod(X*sqrt(w))+T0b)
      m <- v%*%(T0b%*%beta0+t(sqrt(w)*X)%*%(sqrt(w)*(z-Phi3-Phi4*t)))
      beta <- c(spam::rmvnorm(1,m,v))

      # Update phi3
      priorprec <- as.numeric(1/(Sigmaphi[1,1]-Sigmaphi[1,-1]%*%solve(Sigmaphi[-1,-1])%*%Sigmaphi[-1,1]))*Q # Prior Prec of phi3|phi1,phi2,phi4
      priormean <- diag(n)%x%(Sigmaphi[1,-1]%*%solve(Sigmaphi[-1,-1]))%*%c(t(phimat[,-1]))      # Prior mean of phi3|phi1,phi2,phi4
      prec <- priorprec+as.spam(diag(tapply(w,sid,sum),n,n))
      tmp <- tapply(w*(z-X%*%beta-Phi4*t),sid,sum)
      m <- c(priorprec%*%priormean)+tmp
      if(is.positive.definite(prec%>%as.matrix)) phi3 <- spam::rmvnorm.canonical(1, m, prec)[1,]

      # Center
      phi3 <- phi3-mean(phi3)
      Phi3 <- rep(phi3,nis)

      # Update phi4
      priorprec <- as.numeric(1/(Sigmaphi[2,2]-Sigmaphi[2,-2]%*%solve(Sigmaphi[-2,-2])%*%Sigmaphi[-2,2]))*Q # Prior Prec of phi4|phi1,phi2,phi3
      priormean <- diag(n)%x%(Sigmaphi[2,-2]%*%solve(Sigmaphi[-2,-2]))%*%c(t(phimat[,-2]))      # Prior mean of phi4|phi1,phi2,phi3

      prec <- priorprec+as.spam(diag(tapply(w*t^2,sid,sum),n,n))
      tmp <- tapply(t*w*(z-X%*%beta-Phi3),sid,sum)
      m <- c(priorprec%*%priormean)+tmp
      if(is.positive.definite(prec%>%as.matrix)) phi4 <- spam::rmvnorm.canonical(1, m, prec)[1,]

      # Center
      phi4 <- phi4-mean(phi4)
      Phi4 <- rep(phi4,nis)

      # Update Sigma.phi
      phimat <- cbind(phi3,phi4)
      try({
        Sigmaphi <- riwish(2+n-1,diag(2)+t(phimat)%*%Q%*%phimat)
      })

      # Store
      if (i > nburn & i%%nthin==0) {
        j <- (i-nburn)/nthin
        Beta[j,,chain] <- beta
        R[j,chain] <- r
        # R2[j,chain]<-r2
        Sigphi[j,,chain] <- c(Sigmaphi)
        PHI3[j,,chain] <- phi3
        PHI4[j,,chain] <- phi4
        Eta[j,,chain] <- eta
      }

      # if (i%%10==0) print(paste(chain, "/", nchain,"chain | ",round(i/nsim*100,2),"% completed"))
      # if (i%%10==0) print(paste(chain, "/", nchain,"chain | ",round(i/niter*100,2),"% completed |","Test:",conv.test(R[,chain])))

    }
  }

  list_params <- list(Alpha=NULL, Beta=Beta, R=R, Sigphi=Sigphi, PHI3=PHI3, PHI4=PHI4, Eta1=Eta)
  class(list_params) <- "DCMB"
  return(list_params)
}
