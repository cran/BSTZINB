## ----include = FALSE, warning=FALSE, message=FALSE----------------------------
knitr::opts_chunk$set(
  collapse = TRUE, warning = FALSE, message = FALSE,
  fig.width=4, fig.height=4, fig.align='center', 
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages('BSTZINB')

## ----setup--------------------------------------------------------------------
library(BSTZINB)

## -----------------------------------------------------------------------------
data(county.adjacency)
data(USAcities)

library(dplyr)
library(gt)
library(gtsummary)

IAcities <- USAcities %>% filter(state_id=="IA")
countyname <- unique(IAcities$county_name)
A <- get_adj_mat(county.adjacency,countyname,c("IA"))
m <- apply(A,1,sum)

## -----------------------------------------------------------------------------
set.seed(091017)
n <- nrow(A)			 # Number of spatial units
nis <- rep(16,n) 	 # Number of individuals per county; here it's balanced -- 12 per county per year
# Note: may need to lower proposal variance, s, below as n_i increases 
sid <- rep(1:n,nis)
tid <- rep(1:nis[1],n)
N <- length(sid) 		# Total number of observations

## -----------------------------------------------------------------------------
library(spam)
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

## -----------------------------------------------------------------------------
library(boot)
tpop <- USAcities %>% filter(state_id=="IA") %>% group_by(county_fips) %>% summarise(tpop=population) %>% .[,2] %>% unlist %>% as.numeric
x <- rep(log(tpop),nis)
x <- (x-mean(x))/sd(x)
# set.seed(2023); x = rnorm(N)
t <- tid / max(tid)
X <- cbind(1,x)                       # Design matrix, can add additional covariates (e.g., race, age, gender)
Xtilde <- cbind(1,x,sin(t))
p <- ncol(Xtilde)

## -----------------------------------------------------------------------------
# Binomial part
true.alpha <- rep(0.1,p)
true.alpha[1:2] <- c(-0.25,0.25)
# true.alpha <- c(-0.25,0.25,-0.5,0.25)
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

## -----------------------------------------------------------------------------
y <- rep(0,N)                           # Response
y[u==1] <- rnbinom(N1,size=true.r,mu=mu)# If at risk, draw from NB
pzero <- length(y[y==0])/N             # Proportion of zeros
# pzero
# pstruct0/pzero
# table(y)
# mean(mu)
# quantile(tapply(u,id,sum))         # Number of at-risk observations per county *for all 5 years*
# quantile(y)


this_dat <- data.frame(sid,tid,y,X)
head(this_dat)
tbl_summary(this_dat) %>% modify_header(label = "**Coefficients**")

## -----------------------------------------------------------------------------
USDmapCount(state.sel = c("IA"),
            dat = this_dat,
            scol  = 1,
            tcol  = 2,
            cname = countyname,
            uplim=150)

## -----------------------------------------------------------------------------
# blue histogram
tmp <- table(y)/N*100  # convert to %s (divide by N multiply by 100)
oldpar <- par(no.readonly = TRUE)
par(mar=c(3,3,1,1))
barplot(tmp, ylab="Percent",xlab="Count",col="lightblue")

# sphagetti plot
library(CorrMixed)
# Plot individual profiles + mean
Spaghetti.Plot(Dataset = this_dat, Outcome = y, Id=sid, Time = tid)

par(oldpar)

# Moran's I
library(ape)
library(reshape)
this_dat2 <- cast(this_dat,sid~tid,sum,value="y")
Moran.I(rowMeans(this_dat2),A)
this_dat2 %>% apply(2,function(w) Moran.I(w,A)) %>% lapply(function(list) list$p.value) %>% unlist
Moran.I(rowMeans(this_dat2),A)

## -----------------------------------------------------------------------------
USDmapCount(state.sel = c("IA"),
            dat = this_dat,
            scol  = 1,
            tcol  = 2, tsel = 3,
            cname = countyname) ## Timepoint 3

USDmapCount(state.sel = c("IA"),
            dat = this_dat,
            scol  = 1,
            tcol  = 2, tsel = 6,
            cname = countyname) ## Timepoint 6

USDmapCount(state.sel = c("IA"),
            dat = this_dat,
            scol  = 1,
            tcol  = 2, tsel = 9,
            cname = countyname) ## Timepoint 9

USDmapCount(state.sel = c("IA"),
            dat = this_dat,
            scol  = 1,
            tcol  = 2, tsel = 12,
            cname = countyname) ## Timepoint 12

## -----------------------------------------------------------------------------
res0 <- BNB(y, X, A, nchain=2, niter=50, nburn=10)
glimpse(res0)

res1 <- BZINB(y, X, A, nchain=2, niter=50, nburn=10)
glimpse(res1)


apply(res0$Beta,2,mean)
apply(res1$Beta,2,mean)

## -----------------------------------------------------------------------------
res2 <- BSTNB(y, X, A, nchain=2, niter=50, nburn=10)
glimpse(res2)

## ----fig.width=6,fig.height=3,fig.align='center'------------------------------
res3 <- BSTZINB(y, X, A, LinearT=TRUE, nchain=2, niter=50, nburn=10, nthin=1)
glimpse(res3)


conv.test(res3$Alpha) ## Convergence test for Alpha parameters
conv.test(res3$Beta) ## Convergence test for Beta parameters
conv.test(res3$R) ## Convergence test for R parameters

## Plotting the chains for the parameters
oldpar <- par(no.readonly = TRUE)
par(mfrow=c(1,2));plot(res3$R[,1],main="R");plot(res3$R[,2],main="R")
par(mfrow=c(1,2));plot(res3$Alpha[,1,1],type='o',main="Alpha");plot(res3$Alpha[,1,2],type='o',main="Alpha")
par(oldpar)
## Checking the estimated beta with the true beta
apply(res3$Beta,2,mean)
true.beta

## Computing DIC
compute_ZINB_DIC(y,res3,dim(res3$Beta)[1],2)

## ----fig.width=6,fig.height=3,fig.align='center'------------------------------
res4 <- BSTZINB(y, X, A, LinearT=FALSE, nchain=2, niter=50, nburn=10)
glimpse(res4)

## Convergence tests
conv.test(res4$Alpha) ## Alpha parameters
conv.test(res4$Beta) ## Beta parameters
conv.test(res4$R) ## R parameters
oldpar <- par(no.readonly = TRUE)
par(mfrow=c(1,2));plot(res4$R[,1],main="R");plot(res4$R[,2],main="R")
par(mfrow=c(1,2));plot(res4$Alpha[,1,1],type='o',main="Alpha");plot(res4$Alpha[,1,2],type='o',main="Alpha")
par(oldpar)
## Posterior means for beta and R
apply(res4$Beta,2,mean)
mean(res4$R)

## DIC
compute_ZINB_DIC(y,res4,dim(res4$Beta)[1],2)

## -----------------------------------------------------------------------------
ResultTableSummary(res3)

## -----------------------------------------------------------------------------
ResultTableSummary2(y, X, A, nchain=2, niter=50, nburn=10, nthin=1)

## -----------------------------------------------------------------------------
# Map Visualization : Time-averaged Eta
this_result_dat <- this_dat
this_result_dat$y <- apply(res4$Eta1,2,mean)
this_result_dat %>% head
USDmapCount(state.sel = c("IA"),
            dat = this_result_dat,
            scol  = 1,
            cname = countyname)

## -----------------------------------------------------------------------------
# Map Visualization : Spatio-temporal Random Effect of Eta 
this_result_dat <- this_dat
Phi1 <- rep(apply(res4$PHI1,2,mean),nis)
Phi2 <- rep(apply(res4$PHI2,2,mean),nis)
tt <- this_dat$tid / max(this_dat$tid)
this_result_dat$y <- Phi1+Phi2*tt
USDmapCount(state.sel = c("IA"),
            dat = this_result_dat,
            scol  = 1,
            tcol  = 2,
            cname = countyname)

# Map Visualization : Temporal Fixed Effect of Eta 
this_result_dat <- this_dat
Phi1 <- rep(apply(res4$PHI1,2,mean),nis)
Phi2 <- rep(apply(res4$PHI2,2,mean),nis)
tt <- this_dat$tid / max(this_dat$tid)
this_result_dat$y <- apply(res4$Eta1,2,mean) - Phi1+Phi2*tt
USDmapCount(state.sel = c("IA"),
            dat = this_result_dat,
            scol  = 1,
            tcol  = 2,
            cname = countyname)

## -----------------------------------------------------------------------------
# Map Visualization : Time-averaged probability at risk of disease 
this_result_dat <- this_dat
this_result_dat$y <- inv.logit(apply(res4$Eta1,2,mean))
this_result_dat %>% head
USDmapCount(state.sel = c("IA"),
            dat = this_result_dat,
            scol  = 1,
            tcol  = 2,
            cname = countyname)

# Map Visualization : Time-specified probability at risk of disease 
this_result_dat <- this_dat
this_result_dat$y <- inv.logit(apply(res4$Eta1,2,mean))
this_result_dat %>% head
USDmapCount(state.sel = c("IA"),
            dat = this_result_dat,
            scol  = 1,
            tcol  = 2,
            tsel  = 1,
            cname = countyname)

## -----------------------------------------------------------------------------
# Bar Plot of vn-quantile Probability at Risk of disease
qRankPar(state.set=c("IA"),cname=countyname,bstfit=res3,vn=12,
         cex.title=18, cex.lab=12, cex.legend=12)

# Bar Plot of Top vn number of Probability at Risk of disease
qRankParTop(state.set=c("IA"),cname=countyname,bstfit=res3,vn=12,
            cex.title=18, cex.lab=12, cex.legend=12)

## -----------------------------------------------------------------------------
# Line Plot of Nonlinear Time effect abundance
TimetrendCurve(res4,cname=countyname,vn=3,smooth.mode=FALSE,
               cex.title=18, cex.lab=12, cex.legend=12)

# Line Plot of Nonlinear Time effect abundance (smoothed version)
TimetrendCurve(res3,cname=countyname,vn=5,smooth.mode=TRUE,
               cex.title=18, cex.lab=12, cex.legend=12)

