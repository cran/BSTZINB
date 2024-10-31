#' @title Draw spatial maps of various quantities over regions in the US
#'
#' @description
#' Creates a map of any given quantity (at a selected time or averaged over time) for regions in the US specified by state and county
#'
#' @usage USDmapCount(state.sel,dat,scol,tcol=NULL,tsel=NULL,cname,uplim=NULL)
#'
#' @param state.sel character vector giving the selected states
#' @param dat data frame having named components: y - the necessary quantity (numeric), sid - the region indices, tid - the time indices
#' @param scol column index of the spatial regions
#' @param tcol (optional) column index of the time points
#' @param tsel (optional) selected time point
#' @param cname character vector of county names, must match those in USAcities
#' @param uplim (optional) numeric, upper limit for the given quantity
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
#' @import maps
#' @import dplyr
#' @import viridis
#' @import datasets
#'
#' @return spatial map of the required quantity over the specified region
#'
#' @examples
#' data(simdat)
#' data(county.adjacency)
#' data(USAcities)
#' IAcities <- subset(USAcities,state_id=="IA")
#' countyname <- unique(IAcities$county_name)
#' USDmapCount(state.sel="IA",dat=simdat,scol=1,tcol=2,tsel=150,cname=countyname)
#'
#' @export
#'
USDmapCount <- function(state.sel,dat,scol,tcol=NULL,tsel=NULL,cname,uplim=NULL){

  if(!is.character(state.sel)){stop("state.sel must be character vector")}
  if(!is.numeric(scol)){stop("scol must be numeric")}
  if(!is.null(tcol)){if(!is.numeric(tcol)){stop("tcol must be NULL or numeric")}}
  if(!is.null(tsel)){if(!is.numeric(tsel)){stop("tsel must be NULL or numeric")}}
  if(!is.null(tsel)){if(is.null(tcol)){stop("tcol must be supplied if tsel is not null")}}
  if(!is.character(cname)){stop("cname must be character vector")}
  if(!is.null(uplim)){if(!is.numeric(uplim)){stop("uplim must be NULL or numeric")}}
  if(!is.data.frame(dat)){stop("dat must be a data frame")}
  if(is.null(dat$y)){stop("dat must have an entry named y")}
  if(is.null(dat$sid)){stop("dat must have an entry named sid")}
  if(is.null(dat$tid)){stop("dat must have an entry named tid")}
  if(prod(toupper(state.sel) %in% datasets::state.abb) != 1){stop("State abbreviation does not match")}
  USAcities <- BSTZINB::USAcities
  if(prod(toupper(cname) %in% toupper(USAcities$county_name)) != 1){stop("County names do not match")}

  if(is.null(uplim)){
    ulim <- max(dat$y)
  }else{
    ulim <- uplim
  }
  stateabb <- toupper(state.sel)
  statefull <- tolower(datasets::state.name[which(datasets::state.abb == stateabb)])
  state.map <- map_data("county",region=statefull)
  sidd <- names(dat)[scol]
  tidd <- names(dat)[tcol]
  if(is.null(tsel)){
    cnt <- NULL
    zinb.summary <- dat %>% dplyr::group_by(pick(scol)) %>% dplyr::summarise('cnt'=mean(.data$y))
  }else{
    zinb.summary <- dat %>% dplyr::group_by(pick(scol)) %>% dplyr::filter(get(tidd)==tsel)
    zinb.summary$cnt <- zinb.summary$y
    zinb.summary <- zinb.summary %>% dplyr::select(scol,ncol(zinb.summary))
  }

  zinb.summary$cid <- tolower(cname[unlist(zinb.summary[,sidd])])
  map.df <- dplyr::left_join(state.map, zinb.summary, by = c("subregion"="cid"))
  p <- ggplot2::ggplot(data = map.df, aes(x = .data$long, y = .data$lat, group = .data$group)) +
    ggplot2::geom_polygon(aes(fill = cnt), color="black") +
    ggplot2::scale_fill_gradientn(colours = viridis::viridis(25))+
    ggplot2::theme_bw() +
    ggplot2::labs(fill= "",
                  title = "", x="", y="")
  return(p)

}

#' @title Bar plot for time-averaged log-q estimates over quantile-representative counties (descending order)
#'
#' @description
#' Produce a descending order of bar plot for time-averaged log-q estimates over quantile-representative counties
#'
#' @usage qRankPar(state.set,cname,bstfit,vn=12,
#'                cex.title=18, cex.lab=18, cex.legend=18)
#'
#' @param state.set character vector of set of states on which the the graphics is to be made
#' @param cname character vector of the names of the counties
#' @param bstfit the fitted data for BSTP, BSTNB or BSTZINB
#' @param vn positive integer, number of sample counties to display
#' @param cex.title Positive number to control the size of the text of the main title. Defaults to 18.
#' @param cex.lab Positive number to control the size of the text in the axes labels. Defaults to 18.
#' @param cex.legend Positive number to control the size of the text in the legend. Defaults to 18.
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
#' @import boot
#' @import ggplot2
#' @import datasets
#' @import graphics
#'
#' @return bar graph
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
#' qRankPar(state.set=c("IA"),cname=countyname,bstfit=res3,vn=12,
#'          cex.title=18, cex.lab=12, cex.legend=12)
#' }
#'
#' @export
qRankPar <- function(state.set,cname,bstfit,vn=12,
                    cex.title=18, cex.lab=18, cex.legend=18){

  if(!is.character(state.set)){stop("state.set must be character vector")}
  if(!is.character(cname)){stop("cname must be character vector")}
  if(is.null(bstfit$Eta1)){stop("stfit must have an entry named Eta1")}
  if(is.null(bstfit$PHI1)){stop("stfit must have an entry named PHI1")}
  if(!is.numeric(vn)){stop("vn must be a positive integer")}
  if(vn <= 0){stop("vn must be a positive integer")}
  if(prod(toupper(state.set) %in% datasets::state.abb) != 1){stop("State abbreviation does not match")}
  USAcities <- BSTZINB::USAcities
  if(prod(toupper(cname) %in% toupper(USAcities$county_name)) != 1){stop("County names do not match")}


  ns <- dim(bstfit$PHI1)[2]

  qij.mat <- matrix(boot::inv.logit(apply(bstfit$Eta1,2,mean)),nrow=ns)
  zinb.summary <- data.frame("County"=cname, "m" = rowMeans(qij.mat))
  zinb.summary <- zinb.summary[order(zinb.summary$m),]
  zinb.summary$County <- factor(zinb.summary$County)
  zinb.summary.sample <- zinb.summary[floor(seq(1,ns,length.out=vn)),]

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(1,1),mar=c(3,5,1,1))
  p <- ggplot2::ggplot(zinb.summary.sample,aes(x=reorder(.data$County, .data$m), y=.data$m, fill=.data$County)) + ggplot2::geom_bar(alpha=0.8,stat="identity") +
    ggplot2::xlab("") + ggplot2::ylab("Probability at risk") + ggplot2::ylim(0,1) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "")
  return(p + ggplot2::coord_flip()+ggplot2::theme(axis.title.x = element_text(size=cex.title),
                                                  axis.text.x = element_text(size=cex.lab),
                                                  axis.text.y = element_text(size=cex.lab),
                                                  legend.text = element_text(size=cex.legend)))
}

#' @title Bar plot for time-averaged log-q estimates over top ranking counties (descending order)
#'
#' @description
#' Produce a descending order of bar plot for time-averaged log-q estimates over top ranking counties
#'
#' @usage qRankParTop(state.set,cname,bstfit,vn=12,
#'                    cex.title=18, cex.lab=18, cex.legend=18)
#'
#' @param state.set character vector of set of states on which the the graphics is to be made
#' @param cname character vector of the names of the counties
#' @param bstfit the fitted data for BSTP, BSTNB or BSTZINB
#' @param vn positive integer, number of sample counties to display
#' @param cex.title Positive number to control the size of the text of the main title. Defaults to 18.
#' @param cex.lab Positive number to control the size of the text in the axes labels. Defaults to 18.
#' @param cex.legend Positive number to control the size of the text in the legend. Defaults to 18.
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
#' @import ggplot2
#' @import datasets
#' @import graphics
#'
#' @return bar graph
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
#' qRankParTop(state.set=c("IA"),cname=countyname,bstfit=res3,vn=12,
#'              cex.title=18, cex.lab=12, cex.legend=12)
#' }
#'
#' @export
qRankParTop <- function(state.set,cname,bstfit,vn=12,
                       cex.title=18, cex.lab=18, cex.legend=18){

  if(!is.character(state.set)){stop("state.set must be character vector")}
  if(!is.character(cname)){stop("cname must be character vector")}
  if(is.null(bstfit$Eta1)){stop("stfit must have an entry named Eta1")}
  if(is.null(bstfit$PHI1)){stop("stfit must have an entry named PHI1")}
  if(!is.numeric(vn)){stop("vn must be a positive integer")}
  if(vn <= 0){stop("vn must be a positive integer")}
  if(prod(toupper(state.set) %in% datasets::state.abb) != 1){stop("State abbreviation does not match")}
  USAcities <- BSTZINB::USAcities
  if(prod(toupper(cname) %in% toupper(USAcities$county_name)) != 1){stop("County names do not match")}

  ns <- dim(bstfit$PHI1)[2]

  qij.mat <- matrix(boot::inv.logit(apply(bstfit$Eta1,2,mean)),nrow=ns)
  zinb.summary <- data.frame("County"=cname, "m" = rowMeans(qij.mat))
  zinb.summary <- zinb.summary[order(zinb.summary$m,decreasing = TRUE),]
  zinb.summary$County <- factor(zinb.summary$County)
  zinb.summary.sample <- zinb.summary[c(1:vn),]

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(1,1),mar=c(3,5,1,1))
  p <- ggplot2::ggplot(zinb.summary.sample,aes(x=reorder(.data$County, .data$m), y=.data$m, fill=.data$County)) +
    ggplot2::geom_bar(alpha=0.8,stat="identity") +
    ggplot2::xlab("") +
    ggplot2::ylab("Probability at risk") +
    ggplot2::ylim(0,1) + ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "")

  return(p + ggplot2::coord_flip()+ggplot2::theme(axis.title.x = element_text(size=cex.title),
                                                  axis.text.x = element_text(size=cex.lab),
                                                  axis.text.y = element_text(size=cex.lab),
                                                  legend.text = element_text(size=cex.legend)))
}

#' @title Time-trend curve over the study time domain for counties in the US
#'
#' @description
#' Produce a time-trend curve over the study time domain for counties in the US
#'
#' @usage TimetrendCurve(bstfit,cname,vn=5,smooth.mode=TRUE,
#'                      cex.title=18, cex.lab=18, cex.legend=18)
#'
#' @param bstfit fitted object from BSTP, BSTNB or BSTZINB
#' @param cname character vector of county names to use
#' @param vn positive integer, number of sample counties to use
#' @param smooth.mode logical, should splines be fitted to make it smooth
#' @param cex.title Positive number to control the size of the text of the main title. Defaults to 18.
#' @param cex.lab Positive number to control the size of the text in the axes labels. Defaults to 18.
#' @param cex.legend Positive number to control the size of the text in the legend. Defaults to 18.
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
#' @import datasets
#' @import ggplot2
#' @import dplyr
#' @import splines
#' @importFrom reshape melt
#'
#' @return time-trend curves
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
#' TimetrendCurve(res3,cname=countyname,vn=5,smooth.mode=TRUE,cex.title=18, cex.lab=12, cex.legend=12)
#' }
#'
#' @export
TimetrendCurve <- function(bstfit,cname,vn=5,smooth.mode=TRUE,
                          cex.title=18, cex.lab=18, cex.legend=18){

  if(!is.character(cname)){stop("cname must be character vector")}
  if(is.null(bstfit$Eta1)){stop("stfit must have an entry named Eta1")}
  if(is.null(bstfit$PHI1)){stop("stfit must have an entry named PHI1")}
  if(!is.numeric(vn)){stop("vn must be a positive integer")}
  if(vn <= 0){stop("vn must be a positive integer")}
  if(!is.logical(smooth.mode)){stop("smooth.mode must be TRUE/FALSE")}
  USAcities <- BSTZINB::USAcities
  if(prod(toupper(cname) %in% toupper(USAcities$county_name)) != 1){stop("County names do not match")}

  ns <- dim(bstfit$PHI1)[2]
  nt <- dim(bstfit$Eta1)[2]/ns

  time <- c(1:nt)
  df <- data.frame(matrix(apply(bstfit$Eta1,2,mean),nrow=ns)); rownames(df) <- factor(cname)
  df2 <- data.frame(time,t(df[seq(1,ns,length.out=vn),]))
  if(smooth.mode){
    time <- spline(df2[,1],df2[,2])$x
    df2 <- data.frame(time,apply(df2[,-1],2,function(w) {spline(df2[,1],w)$y}))
  }
  dd <- reshape::melt(df2,c("time"))

  spd.plot <- ggplot2::ggplot(dd, aes(x=.data$time,y=.data$value)) +
    ggplot2::geom_line(aes(colour = .data$variable, group = .data$variable),linewidth=1.2)+
    ggplot2::geom_line(aes(.data$time,0),color="red",linewidth=1.2,lty=2) + ggplot2::ylab("") + ggplot2::xlab("") +
    ggplot2::theme(legend.position = "right",
                   legend.title = element_blank(),
                   legend.text = element_text(size=cex.legend),
                   axis.text.x = element_text(size=cex.lab),
                   axis.text.y = element_text(size=cex.lab),
                   axis.title.x = element_text(size=cex.title))

  return(spd.plot)

}
