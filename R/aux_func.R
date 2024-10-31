
#' @title Adjacency matrix for counties of one or many states in the United States
#'
#' @description
#' Creates the adjacency matrix for the supplied counties within the United States using the available neighborhood information
#'
#' @usage get_adj_mat(county.adjacency,Countyvec,Statevec)
#'
#' @param county.adjacency data frame containing the neighborhood information for the counties of the entire US
#' @param Countyvec character vector containing the names of the counties for which the adjacency matrix is to be computed
#' @param Statevec character vector containing the names of the states the supplied counties belong to
#'
#' @return the corresponding adjacency matrix
#'
#' @examples
#' data(county.adjacency)
#' data(USAcities)
#' IAcities <- subset(USAcities,state_id=="IA")
#' countyname <- unique(IAcities$county_name)
#' A <- get_adj_mat(county.adjacency,countyname,c("IA"))
#'
#' @export
get_adj_mat <- function(county.adjacency,Countyvec,Statevec){
  if(!is.data.frame(county.adjacency)){stop("county.adjacency must be a data frame")}
  if(!is.character(Countyvec)){stop("Countyvec must be a character vector")}
  if(!is.character(Statevec)){stop("Statevec must be a character vector")}
  USAcities <- BSTZINB::USAcities
  if(prod(Countyvec %in% USAcities$county_name) != 1){"County names do not match"}
  if(prod(Statevec %in% USAcities$state_id) != 1){"State abbreviations do not match"}

  county.names1 <- tolower(trimws(unlist(lapply(county.adjacency$Countyname, function(x) unlist(strsplit(as.character(x), ","))[1]))))
  county.names2 <- tolower(trimws(unlist(lapply(county.adjacency$neighborname, function(x) unlist(strsplit(as.character(x), ","))[1]))))
  state.abbs1 <- toupper(trimws(unlist(lapply(county.adjacency$Countyname,
                                              function(x) tail(unlist(strsplit(as.character(x), ",")),n=1)))))
  state.abbs2 <- toupper(trimws(unlist(lapply(county.adjacency$neighborname,
                                              function(x) tail(unlist(strsplit(as.character(x), ",")),n=1)))))

  county.names1.full <- unlist(lapply(1:length(county.names1), function(x) paste(county.names1[x], state.abbs1[x], sep="-")))
  county.names2.full <- unlist(lapply(1:length(county.names2), function(x) paste(county.names2[x], state.abbs2[x], sep="-")))

  county.list <- paste(Countyvec,Statevec,sep="-")
  n.county <- length(county.list)

  county.adj.mat <- matrix(0, nrow=n.county, ncol=n.county)
  for (i in 1:dim(county.adjacency)[1]) {
    inx.county <- which(county.list == county.names1.full[i])
    inx.neighbor.county <- which(county.list == county.names2.full[i])
    county.adj.mat[inx.county, inx.neighbor.county] <- 1
    county.adj.mat[inx.neighbor.county, inx.county] <- 1
  }
  return(county.adj.mat)
}


#' @title Summary Table for a fitted object
#'
#' @description
#' Generates a short summary table for a fitted object using BSTP, BSTNB or BSTZINB function
#'
#' @usage ResultTableSummary(bstfit)
#'
#' @param bstfit fitted object using the function BSTP, BSTNB or BSTZINB
#'
#'
#' @import dplyr
#' @import gtsummary
#'
#' @return summary table
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
#' res3 <- BSTZINB(y, X, A, LinearT=TRUE, nchain=3, niter=100, nburn=20, nthin=1)
#' ResultTableSummary(res3)
#' }
#'
#' @export
ResultTableSummary <- function(bstfit){
  if(is.null(bstfit$Beta)){stop("bstfit must have a named component Beta")}

  betamat <- apply(bstfit$Beta,c(1,2),mean)
  if(is.null(bstfit$Alpha)){
    temp <- data.frame(betamat)
    colnames(temp) <- c(paste("b",colnames(betamat),sep="."))
  }else{
    alphamat <- apply(bstfit$Alpha,c(1,2),mean)
    temp <- data.frame(cbind(alphamat,betamat))
    colnames(temp) <- c(paste("a",colnames(alphamat),sep="."),paste("b",colnames(betamat),sep="."))
  }

  temp %>% tbl_summary(statistic = list(all_continuous() ~ "{mean} ({p5}, {p95})", all_categorical() ~ "{n} ({p}%)"),
                       digits = list(all_continuous() ~ c(3,3))) %>%
    modify_header(label = "**Coefficients**",
                  stat_0 = '**Estimates**') %>%
    modify_caption("**BSTZINB Model Coefficients**") %>%
    bold_labels() %>%
    modify_footnote(
      all_stat_cols() ~ "Point estimates (90% credible intervals)")
}


#' @title Generate a summary table of the outputs all different methods given the data
#'
#' @description
#' Fits BSTP, BSTNB and BSTZINB (with linear or non-linear temporal trend) to a given data and summarizes the results in a table
#'
#' @usage ResultTableSummary2(y,X,A,LinearT=FALSE,
#'                            nchain=3,niter=100,nburn=20,nthin=1)
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
#' @import ggplot2
#' @import dplyr
#' @import gtsummary
#' @import BayesLogit
#' @import spam
#' @import MCMCpack
#' @import msm
#' @import splines
#' @import boot
#' @import gt
#' @importFrom matrixcalc is.positive.definite
#'
#'
#' @return summary tables for the different methods
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
#' ResultTableSummary2(y, X, A, LinearT=TRUE, nchain=3, niter=100, nburn=20, nthin=1)
#' }
#' @export
ResultTableSummary2 <- function(y, X, A, LinearT=FALSE, nchain=3, niter=100, nburn=20, nthin=1){

  if(!is.vector(y)){stop("y must be a vector")}
  if(!is.matrix(X)){stop("X must be a matrix")}
  if(!is.matrix(A)){stop("A must be a matrix")}
  if(!is.logical(LinearT)){LinearT <- as.logical(LinearT)}
  if(nchain < 1){stop("nchain must be a positive integer")}
  if(niter < 1){stop("niter must be a positive integer")}
  if(nburn < 0){stop("nburn must be a non-negative integer")}
  if(nthin < 1){stop("nthin must be a positive integer")}
  y <- as.numeric(y)
  if(min(y,na.rm=T)<0){stop("y must be non-negative")}
  if(!is.numeric(X)){stop("X must be numeric")}
  if(!is.numeric(A)){stop("A must be numeric")}

  stfit4 <- BSTZINB(y, X, A, LinearT=LinearT, nchain, niter, nburn, nthin)
  alphamat <- apply(stfit4$Alpha,c(1,2),mean)
  betamat <- apply(stfit4$Beta,c(1,2),mean)
  temp4 <- data.frame(matrix(NA,nrow(alphamat),ncol(alphamat)+ncol(betamat)+2))
  colnames(temp4) <- c("a.t",paste("a",colnames(alphamat),sep="."),"b.t",paste("b",colnames(betamat),sep="."))
  temp4[,paste("a",colnames(alphamat),sep=".")] <- alphamat
  temp4[,paste("b",colnames(betamat),sep=".")]  <- betamat
  DIC4 <- compute_ZINB_DIC(y,stfit4,(niter-nburn)/nthin,nchain)
  ind4 <- rep(NA,ncol(temp4)); names(ind4) <-  colnames(temp4)
  ind4[paste("a",colnames(alphamat),sep=".")] <- conv.test(stfit4$Alpha)
  ind4[paste("b",colnames(betamat),sep=".")] <- conv.test(stfit4$Beta)


  stfit3 <- BSTNB(y, X, A, nchain, niter, nburn, nthin)
  betamat  <- apply(stfit3$Beta,c(1,2),mean)
  temp3 <- data.frame(matrix(NA,nrow(temp4),ncol(temp4)))
  colnames(temp3) <- colnames(temp4)
  temp3[,paste("b",colnames(betamat),sep=".")] <- betamat
  DIC3 <- compute_NB_DIC(y,stfit3,(niter-nburn)/nthin,nchain)
  ind3 <- rep(NA,length(ind4)); names(ind3) <- names(ind4)
  ind3[paste("b",colnames(betamat),sep=".")] <- conv.test(stfit3$Beta)

  stfit2 <- BZINB(y, X, A, nchain, niter, nburn, nthin)
  alphamat <- apply(stfit2$Alpha,c(1,2),mean)
  betamat  <- apply(stfit2$Beta,c(1,2),mean)
  temp2 <- data.frame(matrix(NA,nrow(temp4),ncol(temp4)))
  colnames(temp2) <- colnames(temp4)
  temp2[,paste("a",colnames(alphamat),sep=".")] <- alphamat
  temp2[,paste("b",colnames(betamat),sep=".")]  <- betamat
  DIC2 <- compute_ZINB_DIC(y,stfit2,(niter-nburn)/nthin,nchain)
  ind2 <- rep(NA,length(ind4)); names(ind2) <- names(ind4)
  ind2[paste("a",colnames(alphamat),sep=".")] <- conv.test(stfit2$Alpha)
  ind2[paste("b",colnames(betamat),sep=".")] <- conv.test(stfit2$Beta)

  stfit1 <- BNB(y, X, A, nchain, niter, nburn, nthin)
  betamat <- apply(stfit1$Beta,c(1,2),mean)
  temp1 <- data.frame(matrix(NA,nrow(temp4),ncol(temp4)))
  colnames(temp1) <- colnames(temp4)
  temp1[,paste("b",colnames(betamat),sep=".")]  <- betamat
  DIC1 <- compute_NB_DIC(y,stfit1,(niter-nburn)/nthin,nchain)
  ind1 <- rep(NA,length(ind4)); names(ind1) <- names(ind4)
  ind1[paste("b",colnames(betamat),sep=".")] <- conv.test(stfit1$Beta)

  tabout <- data.frame(rbind(temp4,temp3,temp2,temp1))
  tabout$var <- c(rep("4BSTZINB",nrow(temp4)),
                  rep("3BSTNB",nrow(temp3)),
                  rep("2BZINB",nrow(temp2)),
                  rep("1BNB",nrow(temp1)))

  table <- tabout%>% tbl_summary(by = var,missing = "no",
                                 digits = list(all_continuous() ~ c(3,3)))%>%
    modify_header(label = "**Coefficients**",
                  stat_1 = '**BNB**',
                  stat_2 = '**BZINB**',
                  stat_3 = '**BSTNB**',
                  stat_4 = '**BSTZINB**') %>%
    modify_footnote(
      all_stat_cols() ~ paste("Point estimates (90% credible intervals)",
                              "DIC1 = ",round(DIC1,1),";",
                              "DIC2 = ",round(DIC2,1),";",
                              "DIC3 = ",round(DIC3,1),";",
                              "DIC4 = ",round(DIC4,1)))

  table$table_body$stat_1[table$table_body$stat_1 == "NA (NA, NA)"] <- NA
  table$table_body$stat_2[table$table_body$stat_2 == "NA (NA, NA)"] <- NA
  table$table_body$stat_3[table$table_body$stat_3 == "NA (NA, NA)"] <- NA
  table$table_body$stat_4[table$table_body$stat_4 == "NA (NA, NA)"] <- NA
  table$table_body$stat_1[table$table_body$stat_1 == "0 (NA%)"] <- NA
  table$table_body$stat_2[table$table_body$stat_2 == "0 (NA%)"] <- NA
  table$table_body$stat_3[table$table_body$stat_3 == "0 (NA%)"] <- NA
  table$table_body$stat_4[table$table_body$stat_4 == "0 (NA%)"] <- NA



  table.fin <- table %>% as_gt %>%
    tab_style(
      style = list(cell_text(color = "darkred")),
      locations = cells_body(columns=.data$stat_4,rows=(ind4==FALSE))
    ) %>%
    tab_style(
      style = list(cell_text(color = "darkred")),
      locations = cells_body(columns=.data$stat_3,rows=(ind3==FALSE))
    )%>%
    tab_style(
      style = list(cell_text(color = "darkred")),
      locations = cells_body(columns=.data$stat_2,rows=(ind2==FALSE))
    )%>%
    tab_style(
      style = list(cell_text(color = "darkred")),
      locations = cells_body(columns=.data$stat_1,rows=(ind1==FALSE))
    )

  return(table.fin)

}

#' @title DIC for BSTZINB fitted objects
#'
#' @description
#' Computes DIC for a BSTZINB fitted object
#'
#' @usage compute_ZINB_DIC(y,bstfit,lastit,nchain)
#'
#' @param y vector of counts, must be non-negative, the response used for fitting a BSTZINB model
#' @param bstfit BSTZINB fitted object
#' @param lastit positive integer, size of the chain used to fit BSTZINB
#' @param nchain positive integer, number of chains used to fit BSTZINB
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
#' @import boot
#'
#' @return DIC value
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
#' res3 <- BSTZINB(y, X, A, LinearT=TRUE, nchain=3, niter=100, nburn=20, nthin=1)
#' compute_ZINB_DIC(y,res3,lastit=(100-20)/1,nchain=3)
#' }
#'
#' @export
compute_ZINB_DIC <- function(y,bstfit,lastit,nchain){

  if(is.null(y)){stop("y must be provided")}
  if(!is.vector(y)){stop("y must be a vector")}
  y <- as.numeric(y)
  if(min(y,na.rm=T)<0){stop("y must be non-negative")}
  if(is.null(bstfit$Eta1)){stop("fit must have a named component Eta1")}
  if(is.null(bstfit$Eta2)){stop("fit must have a named component Eta2")}
  if(is.null(bstfit$R)){stop("fit must have a named component R")}
  if(!is.numeric(lastit) | lastit <= 0){stop("lastit must be a positive integer")}
  if(!is.numeric(nchain) | nchain <= 0){stop("nchain must be a positive integer")}

  computeD.avg <- function(y,bstfit){
    eta1.mean <- apply(bstfit$Eta1,2,mean)
    eta2.mean <- apply(bstfit$Eta2,2,mean)
    r.mean    <- mean(bstfit$R)
    pii <- boot::inv.logit(eta1.mean)
    q <- pmax(0.01,pmin(0.99,1/(1+exp(eta2.mean))))
    I <- bstfit$I[dim(bstfit$I)[1],,nchain]
    dNB <- dnbinom(y[I==1],r.mean,q[I==1],log=T)
    comp1 <- log(1-pii[I==0])
    comp2 <- (log(pii[I==1])+dNB)
    return((-2)*sum(comp2))
  }
  computeD.indiv <- function(y,bstfit,iter,chain){
    eta1 <- bstfit$Eta1[iter,,chain]
    eta2 <- bstfit$Eta2[iter,,chain]
    r <- bstfit$R[iter,chain]
    pii <- boot::inv.logit(eta1)
    q <- pmax(0.01,pmin(0.99,1/(1+exp(eta2))))
    I <- bstfit$I[iter,,chain]
    dNB <- dnbinom(y[I==1],r,q[I==1],log=T)
    comp1 <- log(1-pii[I==0])
    comp2 <- (log(pii[I==1])+dNB)
    return((-2)*sum(comp2))
  }
  Dmat <- matrix(0,lastit,nchain)
  for(iter in 1:lastit){
    for(chain in 1:nchain){
      Dmat[iter,chain] <- computeD.indiv(y,bstfit,iter,chain)
    }
  }

  comp1 <- computeD.avg(y,bstfit)
  comp2 <- mean(Dmat)
  DIC <- comp1 + 2*( comp2 - comp1)
  return(DIC)
}

#' @title DIC for BSTNB or BNB fitted objects
#'
#' @description
#' Computes DIC for a BSTNB or BNB fitted object
#'
#' @usage compute_NB_DIC(y,bstfit,lastit,nchain)
#'
#' @param y vector of counts, must be non-negative, the response used for fitting a BSTNB or BSTP model
#' @param bstfit BSTNB or BNB fitted object
#' @param lastit positive integer, size of the chain used to fit BSTZINB
#' @param nchain positive integer, number of chains used to fit BSTZINB
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
#'
#' @return DIC value
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
#' res2 <- BSTNB(y, X, A, nchain=3, niter=100, nburn=20, nthin=1)
#' compute_NB_DIC(y,res2,lastit=(100-20)/1,nchain=3)
#' }
#'
#' @export
compute_NB_DIC <- function(y,bstfit,lastit,nchain){

  if(is.null(y)){stop("y must be provided")}
  if(!is.vector(y)){stop("y must be a vector")}
  y <- as.numeric(y)
  if(min(y,na.rm=T)<0){stop("y must be non-negative")}
  if(is.null(bstfit$Eta1)){stop("fit must have a named component Eta1")}
  if(is.null(bstfit$R)){stop("fit must have a named component R")}
  if(!is.numeric(lastit) | lastit <= 0){stop("lastit must be a positive integer")}
  if(!is.numeric(nchain) | nchain <= 0){stop("nchain must be a positive integer")}

  computeD.avg <- function(y,bstfit){
    eta.mean <- apply(bstfit$Eta1,2,mean)
    r.mean <- mean(bstfit$R)
    q <- pmax(0.01,pmin(0.99,1/(1+exp(eta.mean))))
    dNB <- dnbinom(y,r.mean,q,log=T)
    return((-2)*sum(dNB))
  }
  computeD.indiv <- function(y,bstfit,iter,chain){
    eta <- bstfit$Eta1[iter,,chain]
    r <- bstfit$R[iter,chain]
    q <- pmax(0.01,pmin(0.99,1/(1+exp(eta))))
    dNB <- dnbinom(y,r,q,log=T)
    return((-2)*sum(dNB))
  }
  Dmat <- matrix(0,lastit,nchain)
  for(iter in 1:lastit){
    for(chain in 1:nchain){
      Dmat[iter,chain] <- computeD.indiv(y,bstfit,iter,chain)
    }
  }

  comp1 <- computeD.avg(y,bstfit)
  comp2 <- mean(Dmat)
  DIC <- comp1 + 2*( comp2 - comp1)
  return(DIC)

}

#' @title convergence test for parameters in the fitted objects
#'
#' @description
#' Conducts a test of convergence for a given parameter in the fitted objects using the posterior samples for the said parameter
#'
#' @usage conv.test(params,nchain=3,thshold=1.96)
#'
#' @param params numeric matrix of dimension 2 (iterations x number of parameters, single chain) or 3 (iterations x number of parameters x chain, multiple chains) of posterior samples
#' @param nchain positive integer, number of chains used to fit BSTZINB, BSTNB or BSTP
#' @param thshold positive scalar, the threshold for testing the convergence. Defaults to 1.96
#'
#' @import coda
#'
#' @return logical vector indicating whether convergence was achieved or not
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
#' res3 <- BSTZINB(y, X, A, LinearT=TRUE, nchain=3, niter=100, nburn=20, nthin=1)
#' conv.test(res3$Alpha,nchain=3)
#' }
#'
#' @export
conv.test <- function(params,nchain=3,thshold=1.96){

  if(!is.numeric(nchain) | nchain <= 0){stop("nchain must be a positive integer")}
  if(!is.numeric(thshold) | thshold <= 0){stop("thshold must be a positive number")}

  if(length(dim(params))==3){
    p <- dim(params)[2]
    out <- rep(NA,p)
    for(colid in 1:p){
      mcmc.temp <- mcmc(params[,colid,])
      testout <- geweke.diag(mcmc.temp) %>% unlist
      out[colid] <- as.logical(prod(abs(testout[1:nchain])<thshold))
    }
    return(out)
  }else{
    mcmc.temp <- mcmc(params)
    testout <- geweke.diag(mcmc.temp) %>% unlist
    return(as.logical(prod(abs(testout[1:nchain])<thshold)))
  }

}
