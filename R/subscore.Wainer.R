#' Estimating true subscores using Wainer's augmentation method
#' @description This function estimates subscores using Wainer's augmentation method (Wainer et. al., 2001). 
#' The central idea of this 
#' procedure is that, the estimation of subscores will be improved by shrinking the individual 
#' observed subscores towards some aggregate values (i.e., group mean subscores). The extent of 
#' the shrinkage depends on the closeness of the subscale being estimated with other subscales as 
#' well as reliabilities of all the subscales. Wainer's augmentation is a multivariate version of 
#' Kelly's formula (Kelly, 1947). For details of Wainer's augmentation subscoring method, please refer to Wainer et al. (2001).
#' @param test.data A list that contains datasets of all subtests and the total test, which can be obtained using function 'data.prep'.
#' @return \item{summary}{It contains statistical summary of the augmented subscores (mean, sd, and reliability).}
#' \item{Augmented.subscores}{It contains augmented subscores that are obtained using Wainer's method.} 
#' @import CTT
#' @import stats
#' @examples  
#'        test.data<-data.prep(scored.data,c(3,15,15,20))
#'         
#'         subscore.Wainer(test.data)
#'         
#'         subscore.Wainer(test.data)$summary
#'         subscore.Wainer(test.data)$subscore.augmented
#' @export
#' @references {
#' Wainer, H., Vevea, J., Camacho, F., Reeve, R., Rosa, K., Nelson, L., Swygert, K., & Thissen, D. (2001). 
#' "Augmented scores - "Borrowing strength" to compute scores based on small numbers of items"
#'  In Thissen, D. & Wainer, H. (Eds.), Test scoring (pp.343 - 387). Mahwah, NJ: Lawrence Erlbaum Associates, Inc. 
#' }

subscore.Wainer <- function(test.data) {
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
  summary.list<-list(mean=MEAN, SD=SD, Augmented.reliability=Aug.reliability)
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
