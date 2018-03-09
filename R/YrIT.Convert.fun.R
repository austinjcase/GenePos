#' YrIT.Convert function
#'
#' @description This function converts yellow rust 0-9 infection types to a numeric scale
#' converts raw stripe rust IT values to H and L values
#	 using the standard 0-10 scale for stripe rust seedling values
#' output is a matrix with race names set as row names
#
#' @param x REQUIRED x is the name of thex is a matrix containing the raw IT values to convert to H or L calls,
#' must be read in as a matrix as given by the example files
#' must be a matrix with the race name set as the row names
#' output is a matrix with race names set as row names
#'
#' @param thresh REQUIRED thresh is the threshold IT for calling an line as high or low lines > or equal to thresh will be called high
#'	thresh can be set to one any nummber between 0 and 9
#'
#' @export
#' YrIT.Convert
#'
#'
#'

YrIT.Convert<- function(x, thresh){
  x <-as.matrix(x)
  xc<-apply(x,2, as.numeric)
  rownames(xc)<-rownames(x)
  x<-xc
  thresh <-as.numeric(thresh)
  x[which(x ==10)] <- 10
  x[which(x ==9)] <- 9
  x[which(x ==8)] <- 8
  x[which(x ==7)] <- 7
  x[which(x ==6)] <- 6
  x[which(x ==5)] <- 5
  x[which(x ==4)] <- 4
  x[which(x ==3)] <- 3
  x[which(x ==2)] <- 2
  x[which(x ==1)] <- 1
  x[which(x ==0)] <- 0
  x[which(x =="NA")] <- NA
  x[which(x =="na")] <- NA
  x[which(x =="")] <- NA
  x[which(x ==" ")] <- NA
  xB<-apply(x,2, as.numeric)
  rownames(xB)<-rownames(x)
  x <-as.matrix(xB)
  x[which(x >= thresh)] <- "H"
  x[which(x <  thresh)] <- "L"
  x<-as.matrix(x)
  xb<-x
  xb[which(xb == "H")] <- 1
  xb[which(xb == "L")] <- 0
  xb<-apply(xb,2, as.numeric)
  rownames(xb)<-rownames(x)
  mins<-colMins(xb, na.rm=T)
  maxes<-colMaxs(xb, na.rm=T)
  if (max(mins) == 1) {selects.min<-cbind(colnames(xb),mins )}
  if (max(mins) == 1) {selects.min<-subset(selects.min,mins=="1")}
  if (max(mins) == 1) {selects.min<-selects.min[,-2]}
  if (max(mins) == 1) {warn.msg.min<-c("No avirulent races on",selects.min)}
  if (max(mins) == 1) warning (warn.msg.min)
  if (min(maxes) == 0) {selects.max<-cbind(colnames(xb),maxes )}
  if (min(maxes) == 0) {selects.max<-subset(selects.max,maxes=="0")}
  if (min(maxes) == 0) {selects.max<-selects.max[,-2]}
  if (min(maxes) == 0) {warn.msg.max<- c("No virulent races on", selects.max )}
  if (min(maxes) == 0) warning (warn.msg.max)
  return(x)
}
