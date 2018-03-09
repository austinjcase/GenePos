#' SrIT.Convert function
#'
#' @description This function converts stem rust 0-4 infection types to a numeric scale
#'  SrIT.Convert(x, thresh)
#'	Function SrIT.Convert(x, thresh) converts raw stem rust IT values first to weighted average
#'	data is then converted to H and L values based on the value of thresh
#'	data must conform to the two IT value system as described by Bowden
#'	can be used for both stem rust and leaf rust using the 0-4 stakman scale
#'	output is a matrix with race names set as row names
#
#' @param x REQUIRED x is the name of thex is a matrix containing the raw IT values to convert to H or L calls,
#' must be read in as a matrix as given by the example files
#' must be a matrix with the race name set as the row names
#' output is a matrix with race names set as row names
#'
#' @param thresh REQUIRED thresh is the threshold IT for calling an line as high or low lines > or equal to thresh will be called high
#'	thresh can be set to one of the following to any of the Bowden IT double values ie "44" "3-2" ";;" "0;0;"
#'
#' @export
#' SrIT.Convert
#'
#'

SrIT.Convert<- function(x, thresh){
  x <-as.matrix(x)
  trim <- function (trim_string) gsub("^\\s+|\\s+$", "", trim_string)
  x<-trim(x)
  thresh <-thresh
  library(matrixStats)

  x[which(x =="44")] <- 10
  x[which(x =="3+3+")] <- 9
  x[which(x =="3+3")] <- 8.7
  x[which(x =="3+3-")] <- 8.3
  x[which(x =="3+2+")] <- 8
  x[which(x =="3+2")] <- 7.7
  x[which(x =="3+2-")] <- 7.3
  x[which(x =="3+1+")] <- 7
  x[which(x =="3+1")] <- 6.7
  x[which(x =="3+1-")] <- 6.3
  x[which(x =="3+;")] <- 6
  x[which(x =="33+")] <- 8.3
  x[which(x =="33")] <- 8
  x[which(x =="33-")] <- 7.7
  x[which(x =="32+")] <- 7.3
  x[which(x =="32")] <- 7
  x[which(x =="32-")] <- 6.7
  x[which(x =="31+")] <- 6.3
  x[which(x =="31")] <- 6
  x[which(x =="31-")] <- 5.7
  x[which(x =="3;")] <- 5.3
  x[which(x =="3-3+")] <- 7.7
  x[which(x =="3-3")] <- 7.3
  x[which(x =="3-3-")] <- 7
  x[which(x =="3-2+")] <- 6.7
  x[which(x =="3-2")] <- 6.3
  x[which(x =="3-2-")] <- 6
  x[which(x =="3-1+")] <- 5.7
  x[which(x =="3-1")] <- 5.3
  x[which(x =="3-1-")] <- 5
  x[which(x =="3-;")] <- 4.7
  x[which(x =="2+3+")] <- 7
  x[which(x =="2+3")] <- 6.7
  x[which(x =="2+3-")] <- 6.3
  x[which(x =="2+2+")] <- 6
  x[which(x =="2+2")] <- 5.7
  x[which(x =="2+2-")] <- 5.3
  x[which(x =="2+1+")] <- 5
  x[which(x =="2+1")] <- 4.7
  x[which(x =="2+1-")] <- 4.3
  x[which(x =="2+;")] <- 4
  x[which(x =="23+")] <- 6.3
  x[which(x =="23")] <- 6
  x[which(x =="23-")] <- 5.7
  x[which(x =="22+")] <- 5.3
  x[which(x =="22")] <- 5
  x[which(x =="22-")] <- 4.7
  x[which(x =="21+")] <- 4.3
  x[which(x =="21")] <- 4
  x[which(x =="21-")] <- 3.7
  x[which(x =="2;")] <- 3.3
  x[which(x =="2-3+")] <- 5.7
  x[which(x =="2-3")] <- 5.3
  x[which(x =="2-3-")] <- 5
  x[which(x =="2-2+")] <- 4.7
  x[which(x =="2-2")] <- 4.3
  x[which(x =="2-2-")] <- 4
  x[which(x =="2-1+")] <- 3.7
  x[which(x =="2-1")] <- 3.3
  x[which(x =="2-1-")] <- 3
  x[which(x =="2-;")] <- 2.7
  x[which(x =="1+3+")] <- 5
  x[which(x =="1+3")] <- 4.7
  x[which(x =="1+3-")] <- 4.3
  x[which(x =="1+2+")] <- 4
  x[which(x =="1+2")] <- 3.7
  x[which(x =="1+2-")] <- 3.3
  x[which(x =="1+1+")] <- 3
  x[which(x =="1+1")] <- 2.7
  x[which(x =="1+1-")] <- 2.3
  x[which(x =="1+;")] <- 2
  x[which(x =="13+")] <- 4.3
  x[which(x =="13")] <- 4
  x[which(x =="13-")] <- 3.7
  x[which(x =="12+")] <- 3.3
  x[which(x =="12")] <- 3
  x[which(x =="12-")] <- 2.7
  x[which(x =="11+")] <- 2.3
  x[which(x =="11")] <- 2
  x[which(x =="11-")] <- 1.7
  x[which(x =="1;")] <- 1.3
  x[which(x =="1;0")] <- 1.3
  x[which(x =="10;")] <- 1.3
  x[which(x =="1-3+")] <- 3.7
  x[which(x =="1-3")] <- 3.3
  x[which(x =="1-3-")] <- 3
  x[which(x =="1-2+")] <- 2.7
  x[which(x =="1-2")] <- 2.3
  x[which(x =="1-2-")] <- 2
  x[which(x =="1-1+")] <- 1.7
  x[which(x =="1-1")] <- 1.3
  x[which(x =="1-1-")] <- 1
  x[which(x =="1-;")] <- 0.7
  x[which(x ==";3+")] <- 3
  x[which(x ==";3")] <- 2.7
  x[which(x ==";3-")] <- 2.3
  x[which(x ==";2+")] <- 2
  x[which(x ==";2")] <- 1.7
  x[which(x ==";2-")] <- 1.3
  x[which(x ==";1+")] <- 1
  x[which(x ==";1")] <- 0.7
  x[which(x =="0;1")] <- 0.7
  x[which(x ==";12")] <- 0.7
  x[which(x ==";1-")] <- 0.3
  x[which(x ==";;")] <- 0
  x[which(x =="0;")] <- 0
  x[which(x =="0;0;")] <- 0
  x[which(x =="00")] <- 0
  x[which(x =="0")] <- 0
  x[which(x =="NA")] <- NA
  x[which(x =="na")] <- NA
  x[which(x =="")] <- NA
  x[which(x ==" ")] <- NA
  xB<-apply(x,2, as.numeric)
  rownames(xB)<-rownames(x)
  x <-as.matrix(xB)
  if(thresh == ";;"){ thresh = 0
  } else if (thresh ==  "0;" ){ thresh = 0
  } else if (thresh ==  "00" ){ thresh = 0
  } else if (thresh ==  "0;0;" ){ thresh = 0
  } else if (thresh ==  ";1-" ){ thresh = 0.3
  } else if (thresh ==  "1-;" ){ thresh = 0.7
  } else if (thresh ==  ";1" ){ thresh = 0.7
  } else if (thresh ==  "1-1-" ){ thresh = 1
  } else if (thresh ==  ";1+" ){ thresh = 1
  } else if (thresh ==  "1-1" ){ thresh = 1.3
  } else if (thresh ==  ";2-" ){ thresh = 1.3
  } else if (thresh ==  "1;" ){ thresh = 1.3
  } else if (thresh ==  "1-1+" ){ thresh = 1.7
  } else if (thresh ==  "11-" ){ thresh = 1.7
  } else if (thresh ==  ";2" ){ thresh = 1.7
  } else if (thresh ==  "1-2-" ){ thresh = 2
  } else if (thresh ==  "1+;" ){ thresh = 2
  } else if (thresh ==  ";2+" ){ thresh = 2
  } else if (thresh ==  "11" ){ thresh = 2
  } else if (thresh ==  "1+1-" ){ thresh = 2.3
  } else if (thresh ==  "1-2" ){ thresh = 2.3
  } else if (thresh ==  ";3-" ){ thresh = 2.3
  } else if (thresh ==  "11+" ){ thresh = 2.3
  } else if (thresh ==  "1-2+" ){ thresh = 2.7
  } else if (thresh ==  "1+1" ){ thresh = 2.7
  } else if (thresh ==  "12-" ){ thresh = 2.7
  } else if (thresh ==  "2-;" ){ thresh = 2.7
  } else if (thresh ==  ";3" ){ thresh = 2.7
  } else if (thresh ==  "2-1-" ){ thresh = 3
  } else if (thresh ==  "1-3-" ){ thresh = 3
  } else if (thresh ==  "1+1+" ){ thresh = 3
  } else if (thresh ==  ";3+" ){ thresh = 3
  } else if (thresh ==  "12" ){ thresh = 3
  } else if (thresh ==  "1+2-" ){ thresh = 3.3
  } else if (thresh ==  "2-1" ){ thresh = 3.3
  } else if (thresh ==  "1-3" ){ thresh = 3.3
  } else if (thresh ==  "12+" ){ thresh = 3.3
  } else if (thresh ==  "2;" ){ thresh = 3.3
  } else if (thresh ==  "2-1+" ){ thresh = 3.7
  } else if (thresh ==  "1-3+" ){ thresh = 3.7
  } else if (thresh ==  "1+2" ){ thresh = 3.7
  } else if (thresh ==  "21-" ){ thresh = 3.7
  } else if (thresh ==  "13-" ){ thresh = 3.7
  } else if (thresh ==  "2-2-" ){ thresh = 4
  } else if (thresh ==  "1+2+" ){ thresh = 4
  } else if (thresh ==  "2+;" ){ thresh = 4
  } else if (thresh ==  "21" ){ thresh = 4
  } else if (thresh ==  "13" ){ thresh = 4
  } else if (thresh ==  "2+1-" ){ thresh = 4.3
  } else if (thresh ==  "1+3-" ){ thresh = 4.3
  } else if (thresh ==  "2-2" ){ thresh = 4.3
  } else if (thresh ==  "21+" ){ thresh = 4.3
  } else if (thresh ==  "13+" ){ thresh = 4.3
  } else if (thresh ==  "2-2+" ){ thresh = 4.7
  } else if (thresh ==  "2+1" ){ thresh = 4.7
  } else if (thresh ==  "1+3" ){ thresh = 4.7
  } else if (thresh ==  "22-" ){ thresh = 4.7
  } else if (thresh ==  "3-;" ){ thresh = 4.7
  } else if (thresh ==  "3-1-" ){ thresh = 5
  } else if (thresh ==  "2-3-" ){ thresh = 5
  } else if (thresh ==  "2+1+" ){ thresh = 5
  } else if (thresh ==  "1+3+" ){ thresh = 5
  } else if (thresh ==  "22" ){ thresh = 5
  } else if (thresh ==  "2+2-" ){ thresh = 5.3
  } else if (thresh ==  "3-1" ){ thresh = 5.3
  } else if (thresh ==  "2-3" ){ thresh = 5.3
  } else if (thresh ==  "22+" ){ thresh = 5.3
  } else if (thresh ==  "3;" ){ thresh = 5.3
  } else if (thresh ==  "3-1+" ){ thresh = 5.7
  } else if (thresh ==  "2-3+" ){ thresh = 5.7
  } else if (thresh ==  "2+2" ){ thresh = 5.7
  } else if (thresh ==  "31-" ){ thresh = 5.7
  } else if (thresh ==  "23-" ){ thresh = 5.7
  } else if (thresh ==  "3-2-" ){ thresh = 6
  } else if (thresh ==  "2+2+" ){ thresh = 6
  } else if (thresh ==  "3+;" ){ thresh = 6
  } else if (thresh ==  "31" ){ thresh = 6
  } else if (thresh ==  "23" ){ thresh = 6
  } else if (thresh ==  "3+1-" ){ thresh = 6.3
  } else if (thresh ==  "2+3-" ){ thresh = 6.3
  } else if (thresh ==  "3-2" ){ thresh = 6.3
  } else if (thresh ==  "31+" ){ thresh = 6.3
  } else if (thresh ==  "23+" ){ thresh = 6.3
  } else if (thresh ==  "3-2+" ){ thresh = 6.7
  } else if (thresh ==  "3+1" ){ thresh = 6.7
  } else if (thresh ==  "2+3" ){ thresh = 6.7
  } else if (thresh ==  "32-" ){ thresh = 6.7
  } else if (thresh ==  "3-3-" ){ thresh = 7
  } else if (thresh ==  "3+1+" ){ thresh = 7
  } else if (thresh ==  "2+3+" ){ thresh = 7
  } else if (thresh ==  "32" ){ thresh = 7
  } else if (thresh ==  "3+2-" ){ thresh = 7.3
  } else if (thresh ==  "3-3" ){ thresh = 7.3
  } else if (thresh ==  "32+" ){ thresh = 7.3
  } else if (thresh ==  "3-3+" ){ thresh = 7.7
  } else if (thresh ==  "3+2" ){ thresh = 7.7
  } else if (thresh ==  "33-" ){ thresh = 7.7
  } else if (thresh ==  "3+2+" ){ thresh = 8
  } else if (thresh ==  "33" ){ thresh = 8
  } else if (thresh ==  "3+3-" ){ thresh = 8.3
  } else if (thresh ==  "33+" ){ thresh = 8.3
  } else if (thresh ==  "3+3" ){ thresh = 8.7
  } else if (thresh ==  "3+3+" ){ thresh = 9
  } else if (thresh ==  "44" ){ thresh = 10
  } else { thresh == "NA"
  }
  x[which(x >= thresh)] <- "H"
  x[which(x <  thresh)] <- "L"
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

