#' Estimating true subscores using Yen's OPI
#' @description This function estimates subscores using Yen's Objective Performance Index (OPI; Yen, 1987). 
#' Yen's OPI (Yen, 1987) is a procedure combining Bayesian method and item response theory 
#' (IRT; Embretson & Reise, 2000; Reckase, 1997). 
#' This method pulls an examinee's performance on a certain objective (i.e., subscale) 
#' towards his/her total test performance in order to get a more stable and precise 
#' objective subscore estimate. 
#' @param test.data A list that contains datasets of all subtests and the total test, which can be obtained using function 'data.prep'.
#' @return \item{summary}{It contains statistical summary of OPI (mean & sd).}
#' \item{OPI}{Estimated OPI values} 
#' @import CTT
#' @import stats
#' @import irtoys
#' @examples  
#'         test.data<-data.prep(scored.data,c(3,15,15,20))
#'         
#'         Yen.OPI(test.data)
#'         
#'         Yen.OPI(test.data)$summary
#'         Yen.OPI(test.data)$OPI
#' @export
#' @references {
#' Embretson, S. E., & Reise, S. P. (2013).
#' "Item response theory". Mahwah, NJ: Lawrence Erlbaum Associates, Inc. 
#' }
#' @references {
#' Reckase, M. D. (1997).
#' "The past and future of multidimensional item response theory". 
#' Applied Psychological Measurement, 21(1), 25-36.  
#' }
#' @references {
#' Yen, W. M. (1987, June).
#' "A Bayesian/IRT index of objective performance". 
#' Paper presented at annual meeting of the Psychometric Society, Montreal, Quebec, Canada.
#' }

Yen.OPI<-function(test.data) {

item.prob <- function(a,b,c,theta){item.prob <- 1/(1+exp(-a*(theta-b)))
	return(item.prob)}
p1 <- function(a,c,p){p1 <- a*(1-p)*(p-c)/(1-c)
	return(p1)}

K <- nrow(test.data$total.test)
J <- length(test.data)-1
n <-rep(NA,J)
n.cases<-rep(NA,J)
for (t in 1:J) {
  n[t]<-dim(test.data[[t]])[2] 
  n.cases[t]<-dim(test.data[[t]])[1] 
} 
n.total<-sum(n)

item.par <- est(test.data$total.test, model="2PL", engine="ltm")$est
th.eap <- eap(test.data$total.test, ip=item.par, qu=normal.qu())[,1]
th.eap <- matrix(th.eap,K,n.total,byrow=F)
a.par <- item.par[,1]
b.par <- item.par[,2]
c.par <- item.par[,3]
a <- matrix(a.par,K,n.total,byrow=T)
b <- matrix(b.par,K,n.total,byrow=T)
c <- matrix(c.par,K,n.total,byrow=T)

exp.item.prob <- item.prob(a=a, b=b, c=c, th.eap)	# prob on each item
exp.p.sub <- vector("list", J)
for (j in 1:J){
	if (j == 1){exp.p.sub[[1]] <- apply(exp.item.prob[,1:n[1]],1,sum)/n[1]}
	else{exp.p.sub[[j]] <- apply(exp.item.prob[,(sum(n[1:(j-1)])+1):sum(n[1:j])],1,sum)/n[j]}}

obs.sub <- vector("list", J)
for (j in 1:J){obs.sub[[j]] <- apply(test.data[[j]],1,sum,na.rm=T)}
Q_func <- function(exp.item.prob, obs.sub, n){
	Q <- matrix(NA, K, J)
	for (j in 1:J){
		Q[,j] <- n[j]*(obs.sub[[j]]/n[j] - exp.p.sub[[j]])^2/(exp.p.sub[[j]]*(1-exp.p.sub[[j]]))}
	Q <- apply(Q,1,sum,na.rm=T)
	return(Q)}
  Q <- Q_func(exp.item.prob, obs.sub, n)
  Q.crt <- rep(qchisq(.90, df=J),K)		
  p1 <- p1(a,c,exp.item.prob)
  Info.item <- p1^2/(exp.item.prob*(1-exp.item.prob))
  Info <- apply(Info.item,1,sum,na.rm=T)

sigma2 <- vector("list", J)		
for (j in 1:J){
	if (j == 1){ sigma2[[1]] <- (apply(p1[,1:n[1]],1,sum)/n[1])^2/Info }
	else{ sigma2[[j]] <- (apply(p1[,(sum(n[1:(j-1)])+1):sum(n[1:j])],1,sum)/n[j])^2/Info }
}
n.star <- vector("list", J)		
for (j in 1:J){
	n.star[[j]] <- exp.p.sub[[j]]*(1-exp.p.sub[[j]])/sigma2[[j]] - 1}

OPI <- matrix(NA, K, J)
for (j in 1:J){
	p <- exp.p.sub[[j]] * n.star[[j]] + obs.sub[[j]]
	q <- (1-exp.p.sub[[j]])*n.star[[j]] + n[j] - obs.sub[[j]]

	for (k in 1:K) {			# the k-th examinee
		if (Q[k] <= Q.crt[k]) {
			OPI[k,j] = p[k]/(p[k]+q[k])	
		} else {
   		OPI[k,j] = obs.sub[[j]][k]/n[j]	}}}

colnames(OPI)<-paste("OPI.",1:J,sep="")
SD<-apply(OPI, 2, sd,na.rm=T)
MEAN<-colMeans(OPI,na.rm=T)
summary.list<-list(mean=MEAN, SD=SD)
summary<-do.call(cbind,summary.list)

mylist.names <- c(paste ('Original.Subscore.',rep(1:J),sep=''),'Total.Score')
subscore.list <- as.list(rep(NA, length(mylist.names)))
names(subscore.list) <- mylist.names
for (t in 1 : (length(test.data)))  {
  subscore.list[[t]]<- rowSums(test.data[[t]],na.rm = T)
}  
subscore.original<-do.call(cbind, subscore.list) 
return(list(summary=summary,subscore.original=subscore.original, OPI=OPI))
}
