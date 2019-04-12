#' Computing Yen's Q3 statistic for unidimensional Rasch, 1-, 2-, and 3-PL logistic IRT models 
#' @description This function calculates Yen's Q3 statistics (Yen, 1984; 1993) for unidimensional Rasch, 1-, 2-, 
#' and 3-PL logistic IRT models to assess the local independence assumption.  
#' @param scored.data Item response data with rows as individuals and columns as items.
#' @param IRT.model IRT model ('Rasch', '1pl', '2pl', or '3pl') to be used.The default option is 2pl.
#' @return \item{Q3}{A matrix of Q3 statistics}
#' @return \item{Q3.weighted}{A matrix of Q3 statistics 
#' as obtained by weighting the residual values to reflect 
#' the number of examinees with each response patten.}
#' @import ltm
#' @import stats
#' @import boot
#' @examples  
#'         Yen.Q3(scored.data,IRT.model="2pl")
#'         
#'         Yen.Q3(scored.data)$Q3
#'         Yen.Q3(scored.data)$Q3.weighted
#' @export
#' @references {
#' Yen, W. M. (1984).
#' "Effects of local item dependence on the fit and equating performance of the three-parameter logistic model."
#' Applied Psychological Measurement, 8(2), 125-145. 
#' } 
#'   
#' @references {
#' Yen, W. M. (1993).
#' "Scaling performance assessments: Strategies for managing local item dependence. " 
#' ournal of educational measurement, 30(3), 187-213. 
#' }

Yen.Q3<-function(scored.data,IRT.model="2pl") {
  n<-dim(scored.data)[1]
  n.items<-dim(scored.data)[2]

  if (IRT.model=="Rasch") {
    fit.rasch<-rasch(scored.data, constraint=cbind(ncol(scored.data)+1,1))
    b.s<-coef(fit.rasch)[,1]
    thetas<-ltm::factor.scores(fit.rasch,method = "EAP")$score.dat$z1
    n.thetas<-length(thetas)
    prob.s.2<-matrix(NA,n.thetas,n.items)
    response<-as.matrix(ltm::factor.scores(fit.rasch,method = "EAP")$score.dat[,1:n.items])
    weights<-ltm::factor.scores(fit.rasch,method = "EAP")$score.dat$Obs
    for (j in 1:n.items) {
      for (i in 1:n.thetas) {
        prob.s.2[i,j]<-exp(thetas[i]-b.s[j])/(1+exp(thetas[i]-b.s[j]))
      }}
  }else if (IRT.model=="1pl") {
    fit.1pl<-rasch(scored.data)
    b.s<-coef(fit.1pl)[,1]
    a.s<-coef(fit.1pl)[,2]
    thetas<-ltm::factor.scores(fit.1pl,method = "EAP")$score.dat$z1
    n.thetas<-length(thetas)
    prob.s.2<-matrix(NA,n.thetas,n.items)
    response<-as.matrix(ltm::factor.scores(fit.1pl,method = "EAP")$score.dat[,1:n.items])
    weights<-ltm::factor.scores(fit.1pl,method = "EAP")$score.dat$Obs
        for (j in 1:n.items) {
      for (i in 1:n.thetas) {
        prob.s.2[i,j]<-exp(a.s[j]*(thetas[i]-b.s[j]))/(1+exp(a.s[j]*(thetas[i]-b.s[j])))
      }}
  }else if (IRT.model=="2pl") {
    fit.2pl<-ltm(scored.data~z1)
    b.s<-coef(fit.2pl)[,1]
    a.s<-coef(fit.2pl)[,2]
    thetas<-ltm::factor.scores(fit.2pl,method = "EAP")$score.dat$z1
    n.thetas<-length(thetas)
    prob.s.2<-matrix(NA,n.thetas,n.items)
    response<-as.matrix(ltm::factor.scores(fit.2pl,method = "EAP")$score.dat[,1:n.items])
    weights<-ltm::factor.scores(fit.2pl,method = "EAP")$score.dat$Obs
    for (j in 1:n.items) {
      for (i in 1:n.thetas) {
        prob.s.2[i,j]<-exp(a.s[j]*(thetas[i]-b.s[j]))/(1+exp(a.s[j]*(thetas[i]-b.s[j])))
      }}
  }else if (IRT.model=="3pl") {
    fit.3pl<-tpm(scored.data)
    b.s<-coef(fit.3pl)[,2]
    a.s<-coef(fit.3pl)[,3]
    c.s<-coef(fit.3pl)[,1]
    thetas<-ltm::factor.scores(fit.3pl,method = "EAP")$score.dat$z1
    n.thetas<-length(thetas)
    prob.s.2<-matrix(NA,n.thetas,n.items)
    response<-as.matrix(ltm::factor.scores(fit.3pl,method = "EAP")$score.dat[,1:n.items])
    weights<-ltm::factor.scores(fit.3pl,method = "EAP")$score.dat$Obs
    for (j in 1:n.items) {
      for (i in 1:n.thetas) {
        prob.s.2[i,j]<-c.s[j]+(1-c.s[j])*exp(a.s[j]*(thetas[i]-b.s[j]))/(1+exp(a.s[j]*(thetas[i]-b.s[j])))
      }}
  }
    resid.s<-response-prob.s.2
    Q3.orig<-cor(resid.s)
    Q3.weighted<-matrix(NA,n.items,n.items)
    for (j in 1:n.items) {
      for (k in 1:n.items) {
        Q3.weighted[j,k]<-boot::corr(cbind(resid.s[,j],resid.s[,k]),w=weights)
      }}
    
    
      cat (paste ("There are", sum(Q3.orig[upper.tri(Q3.orig)]>=0.2,na.rm=T), 
                  "Q3 statistics exceed 0.2.",sep=" "))
      cat (paste ("There are",sum(Q3.weighted[upper.tri(Q3.weighted)]>=0.2,na.rm=T), 
                  "Q3 weighted statistics exceed 0.2.",sep=" "))
return (list(Q3=Q3.orig, Q3.weighted=Q3.weighted))}
