#' Computing subscores using Haberman's method based on observed total scores.
#' @description This function estimates true subscores based on observed total scores 
#' using the method introduced by Haberman (2008) <doi:10.3102/1076998607302636>.
#' @param test.data A list that contains item responses of all subtests and 
#' the entire test, which can be obtained using function 'data.prep'.
#' @return \item{summary}{Summary of obtained subscores (e.g., mean, sd).}
#' \item{PRMSE}{PRMSEs of obtained subscores (for Haberman's methods only).}
#' \item{subscore.original}{Original observed subscores and total score.} 
#' \item{subscore.x}{Subscores that are estimated based on the observed total score.}
#' @import CTT
#' @import stats
#' @examples 
#'        test.data<-data.prep(scored.data,c(3,15,15,20), 
#'                             c("Algebra","Geometry","Measurement", "Math"))
#'        
#'        subscore.x(test.data) 
#'        
#'        subscore.x(test.data)$summary
#'        subscore.x(test.data)$PRMSE
#'        subscore.x(test.data)$Correlation
#'        subscore.x(test.data)$Disattenuated.correlation
#'        subscore.x(test.data)$subscore.x
#' @export
#' @references {
#' Haberman, S. J. (2008). 
#' "When can subscores have value?."
#'  Journal of Educational and Behavioral Statistics, 33(2), 204-229.doi:10.3102/1076998607302636
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

  subscore.list <- as.list(rep(NA, n.tests))
  names(subscore.list) <- names(test.data)
  for (t in 1 : (n.tests))  {
    subscore.list[[t]]<- rowSums(test.data[[t]],na.rm = T)
  }  
  
  subscore.original.matrix<-do.call(cbind, subscore.list) 
  corr<-cor(subscore.original.matrix)
  
  for (r in 1:(n.tests)) {
    reliability.alpha[r]<-CTT::itemAnalysis(test.data[[r]],,NA.Delete=T, itemReport=F)$alpha
  } 
  disattenuated.corr<-CTT::disattenuated.cor(corr, reliability.alpha)[-n.tests,-n.tests]
  
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
  
  PRMSE.x<-rep(NA,n.tests)
  r.StXt<-rep(NA,n.tests)
  
  cov.rowsum<-rowSums(CovMat.true[,1:n.subtests],na.rm = TRUE)
  
  for (t in 1:n.subtests) {    
    r.StXt[t]<-cov.rowsum[t]^2/(var.true[t]*var.true[n.tests])
    PRMSE.x[t]<-r.StXt[t]*reliability.alpha[n.tests]
  } 
  mylist.names <- c(paste('Subscore.x.',rep(names(test.data)[-length(test.data)]),sep=''))
  subscore.list.RegOnTot <- as.list(rep(NA, length(mylist.names)))
  names(subscore.list.RegOnTot) <- mylist.names
  
  for (t in 1:n.subtests) { 
    subscore.list.RegOnTot[[t]]<-mean[t]+sqrt(PRMSE.x[t])*(sigma.true[t]/(sigma.obs[n.tests])*(subscore.dataframe[,n.tests]-mean[n.tests]))
  } 
  
  subscore.information.list<-list(Original.reliability=reliability.alpha,PRMSE.x=PRMSE.x)
  subscore.information<-do.call(cbind,subscore.information.list)
  rownames(subscore.information)<-names(test.data)
  
  subscore.original<-do.call(cbind,subscore.list)
  subscore.x<-do.call(cbind,subscore.list.RegOnTot)

  Orig.mean<-mean[1:n.subtests]
  Orig.sd<-SD[1:n.subtests]
  subscore.x.mean<-colMeans(subscore.x,na.rm=T)
  subscore.x.sd<-apply(subscore.x, 2, sd,na.rm=T)

  summary.list<-list(Orig.mean=Orig.mean, Orig.sd=Orig.sd,
                     subscore.x.mean=subscore.x.mean,subscore.x.sd=subscore.x.sd)
  summary<-do.call(cbind,summary.list)
  rownames(summary)<-names(test.data)[1:n.subtests]
  
  if (sum(PRMSE.x[PRMSE.x>1],na.rm=T)>=1) {
    warning ("PRMSE value(s) exceeds 1. The corresponding (augmented) subscore does not have added value.",
             call. = FALSE)
  }
  
  return (list(summary=summary,
               Correlation=corr,
               Disattenuated.correlation=disattenuated.corr, 
               PRMSE=subscore.information, 
               subscore.original=subscore.original,
               subscore.x=subscore.x))  } 
