#' Computing correlation indices for subscores and the total score.
#' @description This function computes Cronback's Alpha and Stratified Alpha. 
#' Disattenuated correlations are also provided.
#' @param test.data A list that contains subscale responses and the total test responses. It 
#' can be obtained using the function 'data.prep'.
#' @return \item{summary}{Summary of obtained subscores (e.g., mean, sd).}
#' \item{correlation}{Correlation indices as indicated above.}
#' @import CTT
#' @import stats
#' @import sirt
#' @examples 
#' # Transfering scored response data to the requried list format
#' test.data<-data.prep(scored.data,c(3,15,15,20),
#'                      c("Algebra","Geometry","Measurement", "Math"))
#'   
#' #Estimate true subscores using Hamerman's method based on observed subscores     
#' subscore.corr(test.data) 
#'        
#' subscore.s(test.data)$summary
#' subscore.s(test.data)$correlation
#' @export
#' @references {
#' Cronbach, L., Schonenman, P., & McKie, D. (1965). 
#' "Alpha coefficients for stratified-parallel tests."
#'  Educational and Psychological Measurement, 25, 291-282. 
#' }

subscore.corr<-function (test.data) {
  
  n.tests<-length(test.data)
  n.subtests<-n.tests-1
  n.items<-rep(NA,n.tests)
  n.cases<-rep(NA,n.tests)
  
  for (t in 1:n.tests) {
    n.items[t]<-dim(test.data[[t]])[2] 
    n.cases[t]<-dim(test.data[[t]])[1] 
  } 
  n.items.total<-n.items[n.tests]
  alpha<-rep(NA, (n.tests))  
  stratefied.alpha<-rep(NA, (n.tests))
  
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
    alpha[r]<-itemAnalysis(test.data[[r]],itemReport=F,NA.Delete=T)$alpha
  } 
  disattenuated.corr<-disattenuated.cor(corr, alpha)[-n.tests,-n.tests]
  Reliabilities<-cbind(alpha, stratefied.alpha)
  rownames(Reliabilities) <- names(test.data)
  return (list(Correlation=corr,
               Disattenuated.correlation=disattenuated.corr, 
               Reliabilities=Reliabilities))}  
