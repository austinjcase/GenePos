#' GenePos function
#' @description  returns the postulated genes and the name of the line Function imports the package "matrixStats"
#
#' @param x REQUIRED is sequential order of the line to test, or the name with "" around it if line.name=TRUE
#' @param lines_HL REQUIRED is the data frame for lines converted to H and L
#' @param isolines_HL REQUIRED is the data frame for isolines converted to H and L
#' @param lines_IT REQUIRED is the data frame for isolines with raw IT values
#' @param isolines_IT REQUIRED is the data frame for lines with raw IT values
#' @param line.name OPTIONAL Defulat is FALSE, if set to TRUE then can specify the name of the line to select
#'
#' @export
#' GenePos


GenePos <- function(x, lines_HL, isolines_HL, isolines_IT, lines_IT, line.name=FALSE){
  y<-lines_HL
  z<-isolines_HL
  w<-isolines_IT
  v<-lines_IT
  y<-as.matrix(y)
  z<-as.matrix(z)
  w<-as.matrix(w)
  wb<-apply(w,2, as.character)
  rownames(wb)<-rownames(w)
  w<-wb
  v<-as.matrix(v)
  vb<-apply(v,2, as.character)
  rownames(vb)<-rownames(v)
  v<-vb
  ifelse(line.name=="FALSE", x<-x, x<-which(colnames(lines_HL)== x))
  if(x > ncol(lines_HL)) warning ("Line is not in data set, value to x wrong")
  library(matrixStats)

  a<-merge(y, z,by=0, all=F) #will do a interetion join for ahmad data
  colnames(a)[1]<-"race"
  b<-a
  b<-as.matrix(b)
  genotypes<-colnames(y)
  b<-b[, c(x+1,(length(genotypes)+2):ncol(b))]
  b[is.na(b)]<-0
  for (i in 1:nrow(b)){
    b[i,2:ncol(b)][which(b[i,1]=="L" & b[i,2:ncol(b)]=="L")]<-1
    b[i,2:ncol(b)][which(b[i,1]=="L" & b[i,2:ncol(b)]=="H")]<-NA
    b[i,2:ncol(b)][which(b[i,1]=="H" & b[i,2:ncol(b)]=="H")]<-NA
    b[i,2:ncol(b)][which(b[i,1]==0 & b[i,2:ncol(b)]=="H")]<-NA
    b[i,2:ncol(b)][which(b[i,1]==0 & b[i,2:ncol(b)]=="L")]<-NA
    b[i,2:ncol(b)][which(b[i,1]=="L" & b[i,2:ncol(b)]==0)]<-NA
    b[i,2:ncol(b)][which(b[i,1]=="H" & b[i,2:ncol(b)]== 0)]<-NA
    b[i,2:ncol(b)][which(b[i,1]==0 & b[i,2:ncol(b)]== 0)]<-NA
    b[i,2:ncol(b)][which(b[i,1]=="H" & b[i,2:ncol(b)]=="L")]<--1
  }
  rownames(b)<-a[,1]
  c<-as.data.frame(b[,-1], stringsAsFactors=F)
  c<-data.matrix(c)
  rownames(c)<-rownames(b)
  colmin<-colMins(c, na.rm=T )
  colmin[which(colmin == "Inf")]<-NA
  rowmin<-rowMins(c, na.rm=T )
  rowmin[which(rowmin == "Inf")]<-NA
  c2<-cbind(rowmin, rownames(c))
  c2<-subset(c2, rowmin >0)
  races<-c2[,2]
  ifelse(	length(races)==0, no.avir<-"TRUE", no.avir<-"FALSE")
  ifelse(length(races)== nrow(isolines_HL), no.vir<-"TRUE", no.vir<-"FALSE")
  c3<-cbind(colmin, colnames(c))
  c4<-1
  c3[is.na(c3)]<-0
  ifelse(max(c3[,1], na.rm=T)<0, c4<-0,  c3<-c3[c3[,1] != -1  ,  ])
  ifelse(length(c3) == 2, genes.one<-"TRUE", genes.one<-"FALSE")
  ifelse(c4>0, ifelse(genes.one=="TRUE", genes<-c3[2], genes<-c3[,2]), genes<-"none")
  ans.gene<-"old"
  name<-genotypes[x]
  line.score<-v[,name]
  ifelse(c4>0, ans<-w[,genes], ans.gene<-"new")
  if(genes.one =="TRUE"){ans<-as.matrix(ans)}
  if(genes.one =="TRUE"){colnames(ans)<-name}
  ifelse(c4>0, boo.select<- rownames(ans) %in% races, ans.gene<-"new")
  ifelse(c4>0,c4b<-cbind(ans, boo.select), ans.gene<-"new")
  ifelse(c4>0,c4b<-subset(c4b, boo.select == "TRUE"), ans.gene<-"new")
  ifelse(c4>0,ans<-c4b, ans.gene<-"new")
  ifelse(ans.gene =="new", ans<-"new", ans<-merge(line.score, ans, by=0 ))
  ifelse(ans.gene =="new", ans<-"new",  colnames(ans)[1:2]<-c("race", name))
  ifelse(genes.one =="TRUE", colnames(ans)[3]<-genes,	ans<-ans)
  ifelse(ans.gene =="new", ans<-"new", row.names(ans)<-ans[,1] )
  ifelse(ans.gene =="new", ans<-"new", ans<-ans[, -1] )
  ifelse(ans.gene =="new", line.score <-as.matrix(line.score), ans<-ans )
  ifelse(ans.gene =="new", races <-as.matrix(races),ans<-ans )
  ifelse(ans.gene =="new",ans.2<-as.matrix(line.score[races,],colnames(junk,)) ,ans<-ans )
  ifelse(ans.gene =="new",colnames(ans.2)<-name, ans<-ans )
  ifelse(ans.gene =="new",ans<-ans.2, ans<-ans )
  ans<-as.data.frame(ans)
  ifelse(ans.gene== "new", ans<-ans, ans<-ans[,-ncol(ans)])
  ifelse(no.avir =="TRUE", ans<-c(name), ans<-ans)
  if( no.avir =="TRUE") warning (c("No avirulent races in data set on ", name))
  if(no.vir =="TRUE") warning (c("No virulent races in data set on ", name))
  if(x > ncol(lines_HL)) ans<-"NA"

  xb<-z
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

  genes<-colnames(ans)
  out<-list("genes"= genes, "races"=ans )
  return(out)
}

