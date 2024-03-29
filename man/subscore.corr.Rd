% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subscore.corr.r
\name{subscore.corr}
\alias{subscore.corr}
\title{Computing correlation indices for subscores and the total score.}
\usage{
subscore.corr(test.data)
}
\arguments{
\item{test.data}{A list that contains item responses of all subtests and the entire test, which can
be obtained using function ’data.prep’.}
}
\value{
\item{summary}{Summary of obtained subscores (e.g., mean, sd).}
\item{correlation}{Correlation indices as indicated above.}
}
\description{
This function computes Cronbach's Alpha and Stratified 
Alpha (Cronbach et al., 1965) <doi: 10.1177/001316446502500201>. 
Disattenuated correlations are also provided.
}
\examples{
# Transferring scored response data to the required list format
test.data<-data.prep(scored.data,c(3,15,15,20),
                     c("Algebra","Geometry","Measurement", "Math"))
  
#Estimate true subscores using Haberman's method based on observed subscores     
subscore.corr(test.data) 
       
subscore.s(test.data)$summary
subscore.s(test.data)$correlation
}
\references{
{
Cronbach, L., Schonenman, P., & McKie, D. (1965). 
"Alpha coefficients for stratified-parallel tests."
 Educational and Psychological Measurement, 25, 291-282.
 doi: 10.1177/001316446502500201. 
}
}
