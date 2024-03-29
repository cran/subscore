#' This function prepares data into a required list format 
#' @description This function generates a list of data sets using the scored original data set, 
#' which can be used as objects in subscore computing functions.
#' @param scored.data Original scored data set with rows as individuals and columns as items.
#' @param subtest.infor A numerical vector. The first number indicates the number of subtests,
#'                       followed by numbers of items on each subscale.
#' @param subtest.names Names of the subscales AND the entire test. The default is NULL. If not provided, names of 
#'                       "subtest.1", "subtest.2",..., will be assigned.                       
#' @return A list that contains item responses of all subtests and the entire test. 
#' The list is then used by other functions (e.g., CTTsub) in the package to obtain subscores. 
#' @examples 
#'         subtest.infor<-c(3,15,15,20) 
#'         subtest.names<-c("Algebra","Geometry","Measurement", "Math")
#'         # This math test consists of 3 subtests, which have 15 algebra 
#'         # items, 15 geometry items, and 20 measurement items.
#'         test.data<-data.prep(scored.data, subtest.infor, subtest.names)
#' @export

data.prep<-function (scored.data, subtest.infor, subtest.names=NULL) {
  n.subtests<-subtest.infor[1]
  
  if (is.null(subtest.names)) {
    test.names<-rep(NA,n.subtests+1)
    test.names[n.subtests+1]<-paste('total.test')
      for (t in 1:n.subtests) { 
        test.names[t]<-paste ('subtest.',t,sep='')    
      }
  } else {
    test.names <- subtest.names}
  
  if ((is.null(subtest.names)==F) & (length(subtest.names)!= (n.subtests+1))) {
    warning ("Name of the total test is needed after the subtest names in subtest.names.",
             call. = FALSE)
  }

  test.data <- as.list(rep(NA, length(test.names)))
  names(test.data) <- test.names
  
  test.data[[1]]<-scored.data[,1:subtest.infor[2]]
  test.data[[n.subtests+1]]<-scored.data
  
  for (t in 2:n.subtests) {  
    test.data[[t]] <-scored.data[,(sum(subtest.infor[2:t],na.rm = TRUE)+1):sum(subtest.infor[2:(t+1)],na.rm = TRUE)]
    } 
  return(test.data)
}  
