#' Computing subscores using Haberman's method based on observed subscores.
#' @description This function estimate true subscores based on observed subscores, 
#' using the method introduced by Haberman (2008).
#' @param test.data A list that contains subscale responses and the total test responses. It 
#' can be obtained using the function 'data.prep'.
#' @return \item{summary}{Summary of obtained subscores (e.g., mean, sd).}
#' \item{PRMSE}{PRMSEs of obtained subscores (for Haberman's methods only).}
#' \item{subscore.original}{Original subscores and total score.}
#' \item{subscore.s}{Subscores that are estimated based on the observed subscore.}
#' @import CTT
#' @import stats
#' @examples 
#' # Transfering scored response data to the requried list format
#' test.data<-data.prep(scored.data,c(3,15,15,20),
#'                      c("Algebra","Geometry","Measurement", "Math"))
#'   
#' #Estimate true subscores using Hamerman's method based on observed subscores     
#' subscore.s(test.data) 
#'        
#' subscore.s(test.data)$summary
#' subscore.s(test.data)$Correlation
#' subscore.s(test.data)$Disattenuated.correlation
#' subscore.s(test.data)$PRMSE
#' subscore.s(test.data)$subscore.s
#' @export
#' @references {
#' Haberman, S. J. (2008). 
#' "When can subscores have value?."
#'  Journal of Educational and Behavioral Statistics, 33(2), 204-229. 
#' }



subscore.s<-function (test.data) {
  
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
  mylist.names <- c(paste('Subscore.s.',rep(names(test.data)[-length(test.data)]),sep=''))
  subscore.list.RegOnSub <- as.list(rep(NA, length(mylist.names)))
  names(subscore.list.RegOnSub) <- mylist.names
  subscore.dataframe<-as.data.frame(subscore.original.matrix)
  for (t in 1: n.subtests) {
    subscore.list.RegOnSub[[t]]<-mean[t]+reliability.alpha[t]*(subscore.dataframe[,t]-mean[t])
  } 
  PRMSE.s<-rep(NA,n.tests)
  PRMSE.s[1:n.subtests]<-reliability.alpha[1:n.subtests]
  
  subscore.information.list<-list(Original.reliability=reliability.alpha, 
                                  PRMSE.s=PRMSE.s)
  subscore.information<-do.call(cbind,subscore.information.list)
  
  rownames.list<-c(paste('Subscore.',rep(names(test.data)[-length(test.data)]),sep=''),names(test.data)[length(test.data)])
  rownames(subscore.information)<-rownames.list
  
  subscore.original<-do.call(cbind,subscore.list)
  subscore.s<-do.call(cbind,subscore.list.RegOnSub)
  
  Orig.mean<-mean[1:n.subtests]
  Orig.sd<-SD[1:n.subtests]
  subscore.s.mean<-colMeans(subscore.s,na.rm=T)
  subscore.s.sd<-apply(subscore.s, 2, sd,na.rm=T)
 
  summary.list<-list(Orig.mean=Orig.mean, Orig.sd=Orig.sd,subscore.s.mean=subscore.s.mean,subscore.s.sd=subscore.s.sd)
  summary<-do.call(cbind,summary.list)
  rownames(summary)<-rownames.list[1:n.subtests]
  
  return (list(summary=summary,
               Correlation=corr,
               Disattenuated.correlation=disattenuated.corr, 
               PRMSE=subscore.information, 
               subscore.original=subscore.original,
               subscore.s=subscore.s))  }  
