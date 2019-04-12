#' This main function estimates true subscores using different methods based on original CTT scores.
#' @description This function estimates true subscores using methods introduced in studies of
#' Haberman (2008) <doi: 10.3102/1076998607302636> and Wainer et al. (2001). 
#' Hypothesis tests (i.e., Olkin' Z,Williams's t, and Hedges-Olkin's Z) are used to determine 
#' whether a subscore or an augmented subscore has added value. Codes for the hypothesis tests are from Sinharay (2019) 
#' <doi: 10.3102/1076998618788862>.
#' @param test.data A list that contains item responses of all subtests and the total test, which can be obtained using function 'data.prep'.
#' @param  method Subscore estimation methods. method="Haberman" (by default) represents the three methods proposed by Harberman (2008)
#' <doi: 10.3102/1076998607302636>. 
#' method="Wainer" represents Wainer's augmented method.         
#' @return \item{summary}{Summary of estimated subscores (e.g., mean, sd).}
#' \item{PRMSE}{(a) PRMSE values of estimated subscores (for Haberman's methods only).(b) Decisions on whether subscores have added 
#' value - added.value.s (or added.value.sx) = 1 means subscore.s (or subscore.sx) has added value, and  added.value.s (or added.value.sx) = 0
#' vice versa.}
#' #' \item{PRMSE.test}{All information in PRMSE plus results of hypopthesis testing based on Sinharay (2019) <doi:10.3102/1076998618788862>.}
#' \item{subscore.original}{Original subscores and total score.}
#' \item{estimated.subscores}{Subscores computed using selected method. Three sets of subscores will be returned if method = "Haberman".}
#' @import CTT
#' @import stats
#' @import sirt
#' @import cocor
#' @examples
#' # Transfering original scored data to a list format
#' # that can be used in other functions.
#' test.data<-data.prep(scored.data,c(3,15,15,20),
#'                      c("Algebra","Geometry","Measurement", "Math"))
#' #----------------------------------------------
#' # Estimating subscores using Haberman's methods       
#' CTTsub(test.data,method="Haberman") # Estimating subscores using Haberman's methods 
#' 
#' # Obtaining original correlation for the three methods
#' CTTsub(test.data,method="Haberman")$Correlation  
#' 
#' # Obtaining disattenuated correlation for the three methods
#' CTTsub(test.data,method="Haberman")$Disattenuated.correlation  
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
#' CTTsub(test.data,method="Haberman")$subscore.s 
#' 
#' # Obtaining subscores that are estimated as a function of the observed total score 
#' CTTsub(test.data,method="Haberman")$subscore.x  
#' 
#' # Obtaining subscores that are estimated as a function of 
#' # both the observed subscores and the observed total score.
#' CTTsub(test.data,method="Haberman")$subscore.sx  
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
#'  Journal of Educational and Behavioral Statistics, 33(2), 204-229. doi:10.3102/1076998607302636. 
#' }
#' @references {
#' Sinharay, S. (2019). 
#' "Added Value of Subscores and Hypothesis Testing."
#'  Journal of Educational and Behavioral Statistics, 44(1), 25-44. doi:10.3102/1076998618788862.
#' }
#' @references {
#' Wainer, H., Vevea, J., Camacho, F., Reeve, R., Rosa, K., Nelson, L., Swygert, K., & Thissen, D. (2001). 
#' "Augmented scores - "Borrowing strength" to compute scores based on small numbers of items."
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
  
  subscore.list <- as.list(rep(NA, n.tests))
  names(subscore.list) <- names(test.data)
  for (t in 1 : (n.tests))  {
    subscore.list[[t]]<- rowSums(test.data[[t]],na.rm = T)
  }  
  
  subscore.original.matrix<-do.call(cbind, subscore.list) 
  corr<-cor(subscore.original.matrix)
  itemstrata<-cbind(colnames(test.data[[n.tests]]),rep(1:(n.tests-1),lengths(test.data)[1:(n.tests-1)]))
  quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  } 
  str.alpha<-quiet(stratified.cronbach.alpha(test.data[[n.tests]], 
                                             itemstrata=itemstrata)$alpha.stratified)
  stratefied.alpha<-c(str.alpha[-1],str.alpha[1])
  
  for (r in 1:(n.tests)) {
    reliability.alpha[r]<-itemAnalysis(test.data[[r]],,NA.Delete=T, itemReport=F)$alpha
  } 
  Reliabilities<-cbind(reliability.alpha, stratefied.alpha)
  disattenuated.corr<-disattenuated.cor(corr, reliability.alpha)[-n.tests,-n.tests]

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
  
  added.value.s<-PRMSE.s>PRMSE.x
  added.value.sx<-(PRMSE.sx-pmax(PRMSE.s,PRMSE.x))>(.1*(1-pmax(PRMSE.s,PRMSE.x)))

  subscore.information.list<-list(Alpha=reliability.alpha, 
                                  PRMSE.s=PRMSE.s, PRMSE.x=PRMSE.x, PRMSE.sx=PRMSE.sx,
                                  added.value.s=added.value.s,added.value.sx=added.value.sx)
  subscore.information<-do.call(cbind,subscore.information.list)
  rownames(subscore.information)<-names(test.data)

  subscore.original<-do.call(cbind,subscore.list)
  subscore.s<-do.call(cbind,subscore.list.RegOnSub)
  subscore.x<-do.call(cbind,subscore.list.RegOnTot)
  subscore.sx<-do.call(cbind,subscore.list.RegOnTotSub)
  
  Orig.mean<-mean[1:n.subtests]
  Orig.sd<-SD[1:n.subtests]
  subscore.s.mean<-colMeans(subscore.s,na.rm=T)
  subscore.x.mean<-colMeans(subscore.x,na.rm=T)
  subscore.sx.mean<-colMeans(subscore.sx,na.rm=T)
  subscore.s.sd<-apply(subscore.s, 2, sd,na.rm=T)
  subscore.x.sd<-apply(subscore.x, 2, sd,na.rm=T)
  subscore.sx.sd<-apply(subscore.sx, 2, sd,na.rm=T)
  
  summary.list<-list(Orig.mean=Orig.mean, Orig.sd=Orig.sd,subscore.s.mean=subscore.s.mean,subscore.s.sd=subscore.s.sd,
                     subscore.x.mean=subscore.x.mean,subscore.x.sd=subscore.x.sd,subscore.sx.mean=subscore.sx.mean,
                     subscore.sx.sd=subscore.sx.sd)
  summary<-do.call(cbind,summary.list)
  rownames(summary)<-names(test.data)[1:n.subtests]
  
  #below code from Sinharay (2019)
  # Compute PRMSEs suggested by Sinharay (2013) from Haberman's PRMSEs
  Olkin.Z<-rep(NA,n.tests)
  Williams.t<-rep(NA,n.tests)
  Hedges.Olkin.Z<-rep(NA,n.tests)
  
  PRs<-PRMSE.s[1:n.subtests]*PRMSE.s[1:n.subtests]
  PRx<-PRMSE.x[1:n.subtests]*PRMSE.s[1:n.subtests]
  PRsx<-PRMSE.sx[1:n.subtests]*PRMSE.s[1:n.subtests]
  rsx<-cor(subscore.original)[1:n.subtests,(n.subtests+1)]
  olk<-rep(0,n.subtests)
  wil<-rep(0,n.subtests)
  n<-n.cases[length(test.data)]
  for (j in 1:n.subtests)
    {olk[j]=cocor.dep.groups.overlap(sqrt(PRs[j]),sqrt(PRx[j]),
                                     rsx[j],n)@olkin1967$statistic
  wil[j]=cocor.dep.groups.overlap(sqrt(PRs[j]),
                                  sqrt(PRx[j]),rsx[j],n)@williams1959$statistic}
  
  compsd<-function(r01,r02,r12,r012,ns,n) {
    a2=2*(r02-r12*r01)/(1-r12*r12)
    a1=-r12*a2
    a3=2*(r12*r01*r01+r12*r02*r02-r01*r02*(1+r12^2))/((1-r12^2)**2)
    V=matrix(0,3,3)
    v11s=(1-r01^2)^2/n
    v22s=(1-r02^2)^2/n
    v33s=(1-r12^2)^2/n
    z=rep(0,4)
   for (j in 1:ns) {
    V[1,2]=(0.5*(2*r12[j]-r01[j]*r02[j])*(1-r12[j]^2-r01[j]^2-r02[j]^2)+r12[j]^3)/n
    V[2,3]=(0.5*(2*r01[j]-r12[j]*r02[j])*(1-r12[j]^2-r01[j]^2-r02[j]^2)+r01[j]^3)/n
    V[1,3]=(0.5*(2*r02[j]-r12[j]*r01[j])*(1-r12[j]^2-r01[j]^2-r02[j]^2)+r02[j]^3)/n
    V[1,1]=v11s[j]
    V[2,2]=v22s[j]
    V[3,3]=v33s[j]
  for (i in 1:2) {for (k in (i+1):3){V[k,i]=V[i,k]}}
    vec=c(a1[j],a2[j],a3[j])
    SD=sqrt(t(vec)%*%V%*%vec)
    z[j]=(r012[j]*r012[j]-r01[j]^2)/SD}
  return(z)}
  
  zs<-compsd(sqrt(PRs),sqrt(PRx),rsx,sqrt(PRsx),n.subtests,n)
  zx<-compsd(sqrt(PRx),sqrt(PRs),rsx,sqrt(PRsx),n.subtests,n)
  hedgesolkin<-ifelse(PRs>PRx,zs,zx)
  
  Olkin.Z[1:n.subtests]<-olk
  Williams.t[1:n.subtests]<-wil
  Hedges.Olkin.Z[1:n.subtests]<-hedgesolkin
  
  subscore.information.list<-list(Alpha=reliability.alpha, 
                                  PRMSE.s=PRMSE.s, PRMSE.x=PRMSE.x, PRMSE.sx=PRMSE.sx,
                                  added.value.s=as.numeric(added.value.s),added.value.sx=as.numeric(added.value.sx),
                                  Olkin.Z=Olkin.Z, Williams.t=Williams.t, 
                                  Hedges.Olkin.Z=Hedges.Olkin.Z)
  subscore.information<-do.call(cbind,subscore.information.list)
  rownames(subscore.information)<-names(test.data)
  
  if (sum(PRMSE.s[PRMSE.s>1],na.rm=T)>=1 | sum(PRMSE.x[PRMSE.x>1],na.rm=T)>=1 | sum(PRMSE.sx[PRMSE.sx>1],na.rm=T)>1) {
    warning ("PRMSE value(s) exceeds 1. The corresponding (augmented) subscore does not have added value.",
             call. = FALSE)
  }
  
  return (list(summary=summary,
               Correlation=corr,
               Disattenuated.correlation=disattenuated.corr, 
               PRMSE=subscore.information[,1:6], 
               PRMSE.test=subscore.information,
               subscore.original=subscore.original,
               subscore.s=subscore.s,
               subscore.x=subscore.x,
               subscore.sx=subscore.sx)) }
  
  else if (method=="Wainer") {
    n.j <- length(test.data)-1
    n.i <- nrow(test.data[[length(test.data)]])
    x <- matrix(NA, n.i, n.j)
    for (j in 1:n.j) {	
      x[, j] <- rowSums(test.data[[j]],na.rm = T)}
    x.bar <- colMeans(x,na.rm = T)
    rho <- c(NA, n.j)
    for (j in 1:n.j){
      rho[j] <- itemAnalysis(test.data[[j]],NA.Delete=T, itemReport=F)$alpha}
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
    colnames(x.true)<-c(paste('Wainer.',rep(names(test.data)[-length(test.data)]),sep=''))
    SD<-apply(x.true, 2, sd,na.rm=T)
    MEAN<-colMeans(x.true,na.rm=T)
    Aug.reliability<-rho.aug
    summary.list<-list(mean=MEAN, SD=SD, Augmented.reliability=Aug.reliability)
    summary<-do.call(cbind,summary.list)
    Augmented.subscores<-x.true
    mylist.names<-names(test.data)
    rownames(summary)<-mylist.names[1:n.j]
    colnames(Augmented.subscores)<-mylist.names[1:n.j]
    subscore.list <- as.list(rep(NA, length(mylist.names)))
    names(subscore.list) <- mylist.names
    for (t in 1 : (length(test.data)))  {
      subscore.list[[t]]<- rowSums(test.data[[t]],na.rm = T)
    }  
    subscore.original<-do.call(cbind, subscore.list) 
    return(list(summary=summary,subscore.original=subscore.original, subscore.augmented=Augmented.subscores))
  }
}  
