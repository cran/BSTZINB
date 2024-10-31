#' @name BNB
#' @title Fit a Bayesian Negative Binomial Model
#'
#' @description
#' Generate posterior samples for the parameters in a Bayesian Negative Binomial Model
#'
#' @usage BNB(y, X, A,
#'            nchain=3, niter=100, nburn=20, nthin=1)
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
#' res0 <- BNB(y, X, A, nchain=2, niter=100, nburn=20, nthin=1)
#' }
#'
#' @export

BNB <- function(y, X, A, nchain=3, niter=100, nburn=20, nthin=1){

  N <- length(y)
  n <- nrow(A)			    # Number of spatial units
  nt <- N/n
  nis <- rep(nt,n) 		# Number of individuals per county; here it's balanced -- 50 per county per year
  # Note: may need to lower proposal variance, s, below as n_i increases
  sid <-rep(1:n,nis)
  tid <-rep(1:nis[1],n)
  N <-length(sid) 		  # Total number of observations
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
  Beta <- array(0,c(lastit,p,nchain))
  colnames(Beta) <- colnames(X)
  R <- matrix(0,lastit,nchain)
  Eta <- array(0,c(lastit,N,nchain))

  for(chain in 1:nchain){

    #########
    # Inits #
    #########
    beta <- alpha <- rnorm(p)
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

      # Update r
      rnew <- rtnorm(1,r,sqrt(s),lower=0)       # Treat r as continuous
      ratio <- sum(dnbinom(y,rnew,q,log=T))-sum(dnbinom(y,r,q,log=T))+
        dtnorm(r,rnew,sqrt(s),0,log=T) - dtnorm(rnew,r,sqrt(s),0,log=T)   # Uniform Prior for R
      # Proposal not symmetric
      if (log(runif(1))<ratio) {
        r <- rnew
        Acc <- Acc+1
      }

      # Update beta
      eta <- X%*%beta
      w <- rpg(N,y+r,eta)                               # Polya weights
      z <- (y-r)/(2*w)
      v <- solve(crossprod(X*sqrt(w))+T0b)
      m <- v%*%(T0b%*%beta0+t(sqrt(w)*X)%*%(sqrt(w)*(z)))
      beta <- c(spam::rmvnorm(1,m,v))

      # Store
      if (i> nburn & i%%nthin==0) {
        j <- (i-nburn)/nthin
        Beta[j,,chain] <- beta
        R[j,chain] <- r
        Eta[j,,chain] <- eta
      }

      # if (i%%10==0) print(paste(chain, "/", nchain,"chain | ",round(i/niter*100,2),"% completed |","Test:",conv.test(R[,chain])))
    }
  }

  list_params <- list(Alpha=NULL, Beta=Beta, R=R, Eta1=Eta)
  return(list_params)
}
