#' Computing subscores using Haberman's method based on both observed total scores and observed subscores.
#' @description This function estimate true subscores based on both observed stotal scores and observed 
#' subscores using the method introduced by Haberman (2008).
#' @param test.data A list that contains subscale responses and the total test responses. It 
#' can be obtained using the function 'data.prep'.
#' @return \item{summary}{Summary of obtained subscores (e.g., mean, sd).}
#' \item{PRMSE}{PRMSEs of obtained subscores (for Haberman's methods only).}
#' \item{subscore.original}{Original observed subscores and total score.} 
#' \item{subscore.sx}{Subscores that are estimated based on both the observed total score and observed subscore.}
#' @import CTT
#' @import stats
#' @examples 
#'        test.data<-data.prep(scored.data,c(3,15,15,20),
#'                             c("Algebra","Geometry","Measurement", "Math"))
#'        
#'        subscore.sx(test.data) 
#'        subscore.s(test.data)$Correlation
#'        subscore.s(test.data)$Disattenuated.correlation
#'        subscore.sx(test.data)$summary
#'        subscore.sx(test.data)$PRMSE
#'        subscore.sx(test.data)$subscore.sx
#' @export
#' @references {
#' Haberman, S. J. (2008). 
#' "When can subscores have value?."
#'  Journal of Educational and Behavioral Statistics, 33(2), 204-229. 
#' }

subscore.sx<-function (test.data) {
  
  n.tests<-length(test.data)
  n.subtests<-n.tests-1
  n.items<-rep(NA,n.tests)
  n.cases<-rep(NA,n.tests)
  
  for (t in 1:n.tests) {
    n.items[t]<-dim(test.data[[t]])[2] 
    n.cases[t]<-dim(test.data[[t]])[1] 
  } 
  n.items.total<-n.items[n.tests]
  reliability.alpha<-rep(NA, (n.tests))  
  
  mylist.names <- names(test.data)
  subscore.list <- as.list(rep(NA, length(mylist.names)))
  names(subscore.list) <- mylist.names
  for (t in 1 : (n.tests))  {
    subscore.list[[t]]<- rowSums(test.data[[t]],na.rm = T)
  }  
  
  subscore.original.matrix<-do.call(cbind, subscore.list) 
  corr<-cor(subscore.original.matrix)
  
  for (r in 1:(n.tests)) {
    reliability.alpha[r]<-itemAnalysis(test.data[[r]],,NA.Delete=T, itemReport=F)$alpha
  } 
  disattenuated.corr<-disattenuated.cor(corr, reliability.alpha)
  if (sum(disattenuated.corr[upper.tri(disattenuated.corr)]>1)>1) {
    warning ("There are disattenuated correlation values exceed 1. PRMSE values should be used with caution 
             and the corresponding (augmented) subscore does not have added value.",
             call. = FALSE)
  }
  sigma.obs<-rep(NA,n.tests)
  for (t in 1:n.tests) {
    sigma.obs[t]<-sd(subscore.list[[t]],na.rm = TRUE)
  }
  
  var.obs<-sigma.obs^2
  CovMat.Obs<-cov(subscore.original.matrix)
  var.true<-var.obs*reliability.alpha
  sigma.true<-sqrt(var.true)
  CovMat.true<-CovMat.Obs
  for (t in 1:n.tests) {
    CovMat.true[t,t]<-var.true[t]
  }
  
  mean<-rep(NA,n.tests)
  SD<-rep(NA,n.tests)
  for (t in 1:n.tests) {
    mean[t]<-mean(subscore.list[[t]],na.rm = TRUE)
    SD[t]<-sd(subscore.list[[t]],na.rm = TRUE)
  }
  subscore.dataframe<-as.data.frame(subscore.original.matrix)
  r.StXt<-rep(NA,n.tests)
  
  cov.rowsum<-rowSums(CovMat.true[,1:n.subtests],na.rm = TRUE)
  
  for (t in 1:n.subtests) {    
    r.StXt[t]<-cov.rowsum[t]^2/(var.true[t]*var.true[n.tests])
  } 
 
  tao<-rep(NA,n.tests)
  beta<-rep(NA,n.tests)
  gamma<-rep(NA,n.tests)
  
  for (t in 1:n.tests) {
    tao[t]<-(sqrt(reliability.alpha[n.tests])*sqrt(r.StXt[t])-corr[t,n.tests]*sqrt(reliability.alpha[t]))/(1-corr[t,n.tests]^2)
    beta[t]<- sqrt(reliability.alpha[t])*(sqrt(reliability.alpha[t])-corr[t,n.tests]*tao[t])
    gamma[t]<-sqrt(reliability.alpha[t])*tao[t]*(sigma.obs[t]/sigma.obs[n.tests])
  } 
  
  PRMSE.sx<-rep(NA, n.tests)
  for (t in 1:n.subtests) { 
    PRMSE.sx[t]<-reliability.alpha[t]+tao[t]^2*(1-corr[t,n.tests]^2)
  } 
  
  mylist.names <- c(paste('Subscore.sx.',rep(names(test.data)[-length(test.data)]),sep=''))
  subscore.list.RegOnTotSub <- as.list(rep(NA, length(mylist.names)))
  names(subscore.list.RegOnTotSub) <- mylist.names
  
  for (t in 1: n.subtests) { 
    subscore.list.RegOnTotSub[[t]]<-mean[t]+beta[t]*(subscore.dataframe[,t]-mean[t])+gamma[t]*(subscore.dataframe[,n.tests]-mean[n.tests])
  } 
  
  subscore.information.list<-list(Original.reliability=reliability.alpha,PRMSE.sx=PRMSE.sx)
  subscore.information<-do.call(cbind,subscore.information.list)
  
  rownames.list<-c(paste('Subscore.',rep(names(test.data)[-length(test.data)]),sep=''),names(test.data)[length(test.data)])
  rownames(subscore.information)<-rownames.list
  
  subscore.original<-do.call(cbind,subscore.list)
  subscore.sx<-do.call(cbind,subscore.list.RegOnTotSub)
  
  Orig.mean<-mean[1:n.subtests]
  Orig.sd<-SD[1:n.subtests]
  subscore.sx.mean<-colMeans(subscore.sx,na.rm=T)
  subscore.sx.sd<-apply(subscore.sx, 2, sd,na.rm=T)
  
  summary.list<-list(Orig.mean=Orig.mean, Orig.sd=Orig.sd,subscore.sx.mean=subscore.sx.mean,
                     subscore.sx.sd=subscore.sx.sd)
  summary<-do.call(cbind,summary.list)
  rownames(summary)<-rownames.list[1:n.subtests]
  
  return (list(summary=summary,
               Correlation=corr,
               Disattenuated.correlation=disattenuated.corr, 
               PRMSE=subscore.information, 
               subscore.original=subscore.original,
               subscore.sx=subscore.sx))  }  
