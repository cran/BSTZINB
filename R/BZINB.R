#' @name BZINB
#' @title Fit a Bayesian Zero Inflated Negative Binomial Model
#'
#' @description
#' Generate posterior samples for the parameters in a Bayesian Zero Inflated Negative Binomial Model
#'
#' @usage BZINB(y,X,A,
#'              nchain=3,niter=100,nburn=20,nthin=1)
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
#' @import BayesLogit
#' @import spam
#' @import MCMCpack
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
#' res1 <- BSTZINB(y, X, A, nchain=2, niter=100, nburn=20, nthin=1)
#' }
#'
#' @export
BZINB <- function(y, X, A, nchain=3, niter=100, nburn=20, nthin=1){

  N <- length(y)
  n <- nrow(A)			    # Number of spatial units
  nt <- N/n
  nis <- rep(nt,n) 		# Number of individuals per county; here it's balanced -- 50 per county per year
  # Note: may need to lower proposal variance, s, below as n_i increases
  sid <- rep(1:n,nis)
  tid <- rep(1:nis[1],n)
  N <- length(sid) 		  # Total number of observations
  p <- ncol(X)
  if(is.null(colnames(X))){colnames(X) <- c("intercept", paste0("X",1:(ncol(X)-1)))}

  ##########
  # Priors #
  ##########
  alpha0 <- beta0 <- rep(0,p)
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
  Beta <- Alpha <- array(0,c(lastit,p,nchain))
  colnames(Beta) <- colnames(Alpha) <- colnames(X)
  R <- matrix(0,lastit,nchain)
  I <- Eta1 <- Eta2 <- array(0,c(lastit,N,nchain))

  for(chain in 1:nchain){

    #########
    # Inits #
    #########
    beta <- alpha <-rnorm(p)
    r <- 1
    Acc <- 0
    y1 <- rep(0,N)             # At risk indicator (this is W in paper)
    y1[y>0] <- 1               # If y>0, then at risk w.p. 1
    N0 <- length(y[y==0])      # Number of observed 0's
    q <- rep(.5,N)             # 1-p=1/(1+exp(X%*%alpha)), used for updating y1

    ########
    # MCMC #
    ########

    for (i in 1:niter){

      # Update alpha
      mu <- X%*%alpha
      w <- rpg(N,1,mu)
      z <- (y1-1/2)/w
      v <- solve(crossprod(sqrt(w)*X)+T0a)
      m <- v%*%(T0a%*%alpha0+t(sqrt(w)*X)%*%(sqrt(w)*z))
      alpha <- c(spam::rmvnorm(1,m,v))

      # Update r
      rnew <- rtnorm(1,r,sqrt(s),lower=0)       # Treat r as continuous
      ratio <- sum(dnbinom(y[y1==1],rnew,q[y1==1],log=T))-sum(dnbinom(y[y1==1],r,q[y1==1],log=T))+
        dtnorm(r,rnew,sqrt(s),0,log=T) - dtnorm(rnew,r,sqrt(s),0,log=T)   # Uniform Prior for R
      # Proposal not symmetric
      if (log(runif(1)) < ratio) {
        r <- rnew
        Acc <- Acc+1
      }

      # Update at-risk indicator y1 (W in paper)
      eta1 <- X%*%alpha
      eta2 <- X%*%beta              # Use all n observations
      pii <- pmax(0.01,pmin(0.99,inv.logit(eta1)))  # at-risk probability
      q <- pmax(0.01,pmin(0.99,1/(1+exp(eta2))))                      # Pr(y=0|y1=1)
      theta <- pii*(q^r)/(pii*(q^r)+1-pii)         # Conditional prob that y1=1 given y=0 -- i.e. Pr(chance zero|observed zero)
      y1[y==0] <- rbinom(N0,1,theta[y==0])      # If y=0, then draw a "chance zero" w.p. theta, otherwise y1=1
      N1 <- sum(y1)
      nis1 <- tapply(y1,sid,sum)

      # Update beta
      eta <- X[y1==1,]%*%beta
      w <- rpg(N1,y[y1==1]+r,eta)                              # Polya weights
      z <- (y[y1==1]-r)/(2*w)                                   # Latent "response"
      v <- solve(crossprod(X[y1==1,]*sqrt(w))+T0b)
      m <- v%*%(T0b%*%beta0+t(sqrt(w)*X[y1==1,])%*%(sqrt(w)*z))
      beta <- c(spam::rmvnorm(1,m,v))

      # Update r2 using Gibbs as in Dadaneh et al and Zhou and Carin #

      # Update latent counts, l
      # l <- rep(0,N1)
      # ytmp <- y[y1==1]
      # for(j in 1:N1) l[j] <- sum(rbinom(ytmp[j],1,round(r/(r+1:ytmp[j]-1),6))) # Could try to avoid loop; rounding avoids numerical stability

      # Update r from conjugate gamma distribution given l and psi
      # eta <- X[y1==1,]%*%beta
      # psi <- exp(eta)/(1+exp(eta))
      # r2 <- rgamma(1,0.01+sum(l),0.01-sum(log(1-psi))) # Gamma(0.01,0.01) prior for r

      # Store
      if (i > nburn & i%%nthin==0) {
        j <- (i-nburn)/nthin
        Alpha[j,,chain] <- alpha
        Beta[j,,chain] <- beta
        R[j,chain] <- r
        # R2[j,chain] <- r2
        Eta1[j,,chain] <- eta1
        Eta2[j,,chain] <- eta2
        I[j,,chain] <- y1
      }

      # if (i%%10==0) print(paste(chain, "/", nchain,"chain | ",round(i/niter*100,2),"% completed |","Test:",conv.test(R[,chain])))
    }
  }

  list_params <- list(Alpha=Alpha, Beta=Beta, R=R,
                     Eta1=Eta1, Eta2=Eta2, I=I)

  return(list_params)
}
