#' Computing subscores using Haberman's method based on observed total scores.
#' @description This function estimates true subscores based on observed total scores 
#' using the method introduced by Haberman (2008).
#' @param test.data A list that contains subscale responses and the total test responses. It 
#' can be obtained using the function 'data.prep'.
#' @return \item{summary}{Summary of obtained subscores (e.g., mean, sd).}
#' \item{PRMSE}{PRMSEs of obtained subscores (for Haberman's methods only).}
#' \item{subscore.original}{Original observed subscores and total score.} 
#' \item{subscore.x}{Subscores that are estimated based on the observed total score.}
#' @import CTT
#' @import stats
#' @examples 
#'        test.data<-data.prep(scored.data,c(3,15,15,20))
#'        
#'        subscore.x(test.data) 
#'        
#'        subscore.x(test.data)$summary
#'        subscore.x(test.data)$PRMSE
#'        subscore.x(test.data)$subscore.x
#' @export
#' @references {
#' Haberman, S. J. (2008). 
#' "When can subscores have value?."
#'  Journal of Educational and Behavioral Statistics, 33(2), 204-229. 
#' }


subscore.x<-function (test.data) {
  
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
  
  mylist.names <- c(paste ('Original.Subscore.',rep(1:n.subtests),sep=''),'Total.Score')
  subscore.list <- as.list(rep(NA, length(mylist.names)))
  names(subscore.list) <- mylist.names
  for (t in 1 : (n.tests))  {
    subscore.list[[t]]<- rowSums(test.data[[t]],na.rm = T)
  }  
  
  subscore.original.matrix<-do.call(cbind, subscore.list) 
  
  for (r in 1:(n.tests)) {
    reliability.alpha[r]<-reliability(test.data[[r]],itemal=TRUE,NA.Delete=T)[[3]]
  } 
  
  sigma.obs<-rep(NA,n.tests)
  for (t in 1:n.tests) {
    sigma.obs[t]<-sd(subscore.list[[t]],na.rm = TRUE)
  }
  
  var.obs<-sigma.obs^2
  corr<-cor(subscore.original.matrix)
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
  
  PRMSE.x<-rep(NA,n.tests)
  r.StXt<-rep(NA,n.tests)
  
  cov.rowsum<-rowSums(CovMat.true[,1:n.subtests],na.rm = TRUE)
  
  for (t in 1:n.subtests) {    
    r.StXt[t]<-cov.rowsum[t]^2/(var.true[t]*var.true[n.tests])
    PRMSE.x[t]<-r.StXt[t]*reliability.alpha[n.tests]
  } 
  mylist.names <- c(paste ('Subscore.x.',rep(1:n.subtests),sep=''))
  subscore.list.RegOnTot <- as.list(rep(NA, length(mylist.names)))
  names(subscore.list.RegOnTot) <- mylist.names
  
  for (t in 1:n.subtests) { 
    subscore.list.RegOnTot[[t]]<-mean[t]+sqrt(PRMSE.x[t])*(sigma.true[t]/(sigma.obs[n.tests])*(subscore.dataframe[,n.tests]-mean[n.tests]))
  } 
  
  subscore.information.list<-list(Original.reliability=reliability.alpha,PRMSE.x=PRMSE.x)
  subscore.information<-do.call(cbind,subscore.information.list)
  
  rownames.list<-c(paste('Subscore.',rep(1:n.subtests),sep=''),'Total.test')
  rownames(subscore.information)<-rownames.list
  
  subscore.original<-do.call(cbind,subscore.list)
  subscore.x<-do.call(cbind,subscore.list.RegOnTot)

  Orig.mean<-mean[1:n.subtests]
  Orig.sd<-SD[1:n.subtests]
  subscore.x.mean<-colMeans(subscore.x,na.rm=T)
  subscore.x.sd<-apply(subscore.x, 2, sd,na.rm=T)

  summary.list<-list(Orig.mean=Orig.mean, Orig.sd=Orig.sd,
                     subscore.x.mean=subscore.x.mean,subscore.x.sd=subscore.x.sd)
  summary<-do.call(cbind,summary.list)
  rownames(summary)<-rownames.list[1:n.subtests]
  
  return (list(summary=summary,
               PRMSE=subscore.information, 
               subscore.original=subscore.original,
               subscore.x=subscore.x))  } 