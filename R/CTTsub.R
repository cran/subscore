#' This main function estimates true subscores using different methods based on original CTT scores.
#' @description This function estimates true subscores using methods introduced in Haberman (2008) and Wainer et al. (2001).
#' @param test.data A list that contains item responses of all subtests and the total test, which can be obtained using function "data.prep".
#' @param  method Subscore estimation methods. method="Haberman" (by default) represents the three methods propose by Harberman (2008). 
#' method="Wainer" represents Wainer's augmented method.         
#' @return \item{summary}{Summary of estimated subscores (e.g., mean, sd).}
#' \item{PRMSE}{PRMSEs of estimated subscores (for Haberman's methods only).}
#' \item{subscore.original}{Original subscores and total score.}
#' \item{estimated.subscores}{Subscores computed using selected methods. Three sets of subscores will be returned if method = "Haberman".}
#' @import CTT
#' @import stats
#' @examples
#' # Transfering original scored data to a list format
#' # that can be used in other functions.
#' test.data<-data.prep(scored.data,c(3,15,15,20))
#' #----------------------------------------------
#' # Estimating subscores using Haberman's methods       
#' CTTsub(test.data,method="Haberman") # Estimating subscores using Haberman's methods 
#' 
#' # Obtaining PRMSEs for the three methods
#' CTTsub(test.data,method="Haberman")$PRMSE  
#' 
#' # Obtaining descriptive statistics summary for estimated subscores  
#' CTTsub(test.data,method="Haberman")$summary 
#' 
#' # Obtaining raw subscores  
#' CTTsub(test.data,method="Haberman")$subscore.original  
#' 
#' # Obtaining subscores that are estimated as a function of the observed subscores 
#' CTTsub(test.data,method="Haberman")$subscore.RegOnSub 
#' 
#' # Obtaining subscores that are estimated as a function of the observed total score 
#' CTTsub(test.data,method="Haberman")$subscore.RegOnTot  
#' 
#' # Obtaining subscores that are estimated as a function of 
#' # both the observed subscores and the observed total score.
#' CTTsub(test.data,method="Haberman")$subscore.RegOnTotSub  
#' 
#' #-------------------------------------------      
#' # Estimating subscores using Wainer's method
#' CTTsub(test.data,method="Wainer") 
#'        
#' # Obtaining descriptive statistics summary for subscores
#' CTTsub(test.data,method="Wainer")$summary   
#' 
#' # Obtaining original subscores
#' CTTsub(test.data,method="Wainer")$subscore.original 
#' 
#' # Obtaining subscores that are estimated using Wainer's augmentation method  
#' CTTsub(test.data,method="Wainer")$subscore.augmented  
#'         
#' @export  
#' @references {
#' Haberman, S. J. (2008). 
#' "When can subscores have value?."
#'  Journal of Educational and Behavioral Statistics, 33(2), 204-229. 
#' }
#' @references {
#' Wainer, H., Vevea, J., Camacho, F., Reeve, R., Rosa, K., Nelson, L., Swygert, K., & Thissen, D. (2001). 
#' "Augmented scores - "Borrowing strength" to compute scores based on small numbers of items"
#'  In Thissen, D. & Wainer, H. (Eds.), Test scoring (pp.343 - 387). Mahwah, NJ: Lawrence Erlbaum Associates, Inc. 
#' }

CTTsub<-function (test.data, method="Haberman") {
  if (method=="Haberman") {
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
  mylist.names <- c(paste ('RegOnSub.Score.',rep(1:n.subtests),sep=''))
  subscore.list.RegOnSub <- as.list(rep(NA, length(mylist.names)))
  names(subscore.list.RegOnSub) <- mylist.names
  subscore.dataframe<-as.data.frame(subscore.original.matrix)
  for (t in 1: n.subtests) {
    subscore.list.RegOnSub[[t]]<-mean[t]+reliability.alpha[t]*(subscore.dataframe[,t]-mean[t])
  } 
  PRMSE.RegOnSub<-rep(NA,n.tests)
  PRMSE.RegOnSub[1:n.subtests]<-reliability.alpha[1:n.subtests]

  PRMSE.RegOnTot<-rep(NA,n.tests)
  r.StXt<-rep(NA,n.tests)
  
  cov.rowsum<-rowSums(CovMat.true[,1:n.subtests],na.rm = TRUE)

  for (t in 1:n.subtests) {    
    r.StXt[t]<-cov.rowsum[t]^2/(var.true[t]*var.true[n.tests])
    PRMSE.RegOnTot[t]<-r.StXt[t]*reliability.alpha[n.tests]
  } 
  mylist.names <- c(paste ('RegOnTot.Score.',rep(1:n.subtests),sep=''))
  subscore.list.RegOnTot <- as.list(rep(NA, length(mylist.names)))
  names(subscore.list.RegOnTot) <- mylist.names
  
  for (t in 1:n.subtests) { 
    subscore.list.RegOnTot[[t]]<-mean[t]+sqrt(PRMSE.RegOnTot[t])*(sigma.true[t]/(sigma.obs[n.tests])*(subscore.dataframe[,n.tests]-mean[n.tests]))
    } 
  
  tao<-rep(NA,n.tests)
  beta<-rep(NA,n.tests)
  gamma<-rep(NA,n.tests)

  for (t in 1:n.tests) {
    tao[t]<-(sqrt(reliability.alpha[n.tests])*sqrt(r.StXt[t])-corr[t,n.tests]*sqrt(reliability.alpha[t]))/(1-corr[t,n.tests]^2)
    beta[t]<- sqrt(reliability.alpha[t])*(sqrt(reliability.alpha[t])-corr[t,n.tests]*tao[t])
    gamma[t]<-sqrt(reliability.alpha[t])*tao[t]*(sigma.obs[t]/sigma.obs[n.tests])
  } 
 
  PRMSE.RegOnTotSub<-rep(NA, n.tests)
  for (t in 1:n.subtests) { 
    PRMSE.RegOnTotSub[t]<-reliability.alpha[t]+tao[t]^2*(1-corr[t,n.tests]^2)
  } 

  mylist.names <- c(paste ('RegOnTotSub.Score.',rep(1:n.subtests),sep=''))
  subscore.list.RegOnTotSub <- as.list(rep(NA, length(mylist.names)))
  names(subscore.list.RegOnTotSub) <- mylist.names

  for (t in 1: n.subtests) { 
    subscore.list.RegOnTotSub[[t]]<-mean[t]+beta[t]*(subscore.dataframe[,t]-mean[t])+gamma[t]*(subscore.dataframe[,n.tests]-mean[n.tests])
  } 
  
  subscore.information.list<-list(Original.reliability=reliability.alpha, 
                                  PRMSE.RegOnSub=PRMSE.RegOnSub, PRMSE.RegOnTot=PRMSE.RegOnTot, PRMSE.RegOnTotSub=PRMSE.RegOnTotSub)
  subscore.information<-do.call(cbind,subscore.information.list)
  
  rownames.list<-c(paste('Subscore.',rep(1:n.subtests),sep=''),'Total.test')
  rownames(subscore.information)<-rownames.list

  subscore.original<-do.call(cbind,subscore.list)
  subscore.RegOnSub<-do.call(cbind,subscore.list.RegOnSub)
  subscore.RegOnTot<-do.call(cbind,subscore.list.RegOnTot)
  subscore.RegOnTotSub<-do.call(cbind,subscore.list.RegOnTotSub)
  
  Orig.mean<-mean[1:n.subtests]
  Orig.sd<-SD[1:n.subtests]
  RegOnSub.mean<-colMeans(subscore.RegOnSub,na.rm=T)
  RegOnTot.mean<-colMeans(subscore.RegOnTot,na.rm=T)
  RegOnTotSub.mean<-colMeans(subscore.RegOnTotSub,na.rm=T)
  RegOnSub.sd<-apply(subscore.RegOnSub, 2, sd,na.rm=T)
  RegOnTot.sd<-apply(subscore.RegOnTot, 2, sd,na.rm=T)
  RegOnTotSub.sd<-apply(subscore.RegOnTotSub, 2, sd,na.rm=T)
  
  summary.list<-list(Orig.mean=Orig.mean, Orig.sd=Orig.sd,RegOnSub.mean=RegOnSub.mean,RegOnSub.sd=RegOnSub.sd,
                     RegOnTot.mean=RegOnTot.mean,RegOnTot.sd=RegOnTot.sd,RegOnTotSub.mean=RegOnTotSub.mean,
                     RegOnTotSub.sd=RegOnTotSub.sd)
  summary<-do.call(cbind,summary.list)
  rownames(summary)<-rownames.list[1:n.subtests]
  
  return (list(summary=summary,
               PRMSE=subscore.information, 
               subscore.original=subscore.original,
               subscore.RegOnSub=subscore.RegOnSub,
               subscore.RegOnTot=subscore.RegOnTot,
               subscore.RegOnTotSub=subscore.RegOnTotSub)) }
  else if (method=="Wainer") {
    n.j <- length(test.data)-1
    n.i <- nrow(test.data$total.test)
    x <- matrix(NA, n.i, n.j)
    for (j in 1:n.j) {	
      x[, j] <- rowSums(test.data[[j]],na.rm = T)}
    x.bar <- colMeans(x,na.rm = T)
    rho <- c(NA, n.j)
    for (j in 1:n.j){
      rho[j] <- reliability(test.data[[j]],NA.Delete=T)$alpha}
    S.obs <- cov(x,use="pairwise.complete.obs")
    S.true <- S.obs
    for (j in 1:n.j) {
      S.true[j,j] = rho[j] * S.obs[j,j]}
    B <- S.true %*% solve(S.obs)
    x.true <- matrix(x.bar, n.j, n.i) + B %*% (t(x) - matrix(x.bar, n.j, n.i))
    x.true <- x.bar + B %*% (t(x) - x.bar)
    x.true <- t(x.true)
    A <- S.true %*% solve(S.obs) %*% S.true %*% solve(S.obs) %*% S.true
    C <- S.true %*% solve(S.obs) %*% S.true
    rho.aug <- c(NA, n.j)
    for (j in 1:n.j){
      rho.aug[j] <- A[j,j]/C[j,j]}
    colnames(x.true)<-paste("subscore.",1:n.j,sep="")
    SD<-apply(x.true, 2, sd,na.rm=T)
    MEAN<-colMeans(x.true,na.rm=T)
    Aug.reliability<-rho.aug
    summary.list<-list(mean=MEAN, sd=SD, original.reliability=rho,augmented.reliability=Aug.reliability)
    summary<-do.call(cbind,summary.list)
    Augmented.subscores<-x.true
    mylist.names <- c(paste ('Original.Subscore.',rep(1:n.j),sep=''),'Total.Score')
    subscore.list <- as.list(rep(NA, length(mylist.names)))
    names(subscore.list) <- mylist.names
    for (t in 1 : (length(test.data)))  {
      subscore.list[[t]]<- rowSums(test.data[[t]],na.rm = T)
    }  
    subscore.original<-do.call(cbind, subscore.list) 
    return(list(summary=summary,subscore.original=subscore.original, subscore.augmented=Augmented.subscores))
  }
} 