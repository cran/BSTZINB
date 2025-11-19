## PRINT METHOD ##

#' @aliases print.DCMB
#' @title Print function for DCMB class of objects
#'
#' @description
#' Prints out the objects of class DCMB
#'
#' @usage \method{print}{DCMB}(x,digits=3,...)
#'
#' @param x object of class DCMB
#' @param digits non-negative integer determining the number of significant digits to print Defaults to 3
#' @param ... additional arguments to pass to the print function
#'
#' @importFrom methods is
#'
#'
#' @return prints out the class object details
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
#' print(res1)
#' }
#'
#' @export print.DCMB
#' @export
print.DCMB <- function(x,digits=3,...){
  if(digits<0){stop("digits must be positive integer")}
  digit <- floor(digits);
  if(!is(x,"DCMB")){stop("x must be an object of class DCMB")}

  nms <- names(x)
  nml <- length(nms)


  for(i in 1:nml){
    cat("Parameter ",nms[i],"\n",sep="")
    dm <- dim(x[[nms[i]]])
    cat("dimension: number of samples=",dm[1]," x number of params=",dm[2]," x number of chains=",dm[3],"\n",sep="")
    print(round(head(x[[nms[i]]]),digits=digit))
    cat("\n")
  }

  invisible()
}

## SUMMARY METHOD ##

#' @aliases summary.DCMB
#' @title Summary function for objects of class DCMB
#'
#' @description
#' Gives out a summary of the posterior samples for parameters of any of the models, outputs of which are contained in a DCMB object
#'
#' @usage \method{summary}{DCMB}(x,...)
#'
#' @param x object of class DCMB
#' @param ... additional parameters to pass on to the function
#'
#' @import dplyr
#' @import gtsummary
#' @importFrom methods is
#'
#' @return returns a table of summary values
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
#' summary(res3)
#' }
#'
#' @export summary.DCMB
#' @export
summary.DCMB <- function(x,...){
  if(!is(x,"DCMB")){stop("x must be an object of class DCMB")}
  vv <- ResultTableSummary(x)
  class(vv) <- "summ.DCMB"
  return(vv)
}


## PRINT SUMMARY METHOD ##

#' @aliases print.summ.DCMB
#' @title Prints out the summary of a DCMB object
#'
#' @description
#' Prints the summary object created by summary function fro DCMB objects
#'
#' @usage \method{print}{summ.DCMB}(x,...)
#'
#' @param x a summary object generated from a DCMB object
#' @param ... additional parameters to pass onto the function
#'
#' @importFrom methods is
#'
#' @return prints the summary of the DCMB object from which the summary object was formed
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
#' print(summary(res3))
#' }
#'
#' @export print.summ.DCMB
#' @export
print.summ.DCMB <- function(x,...){
  if(!is(x,"summ.DCMB")){stop("x must be a summary object generated from a DCMB object")}
  print(x)

  invisible()
}
