% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CTTsub.R
\name{CTTsub}
\alias{CTTsub}
\title{This main function estimates true subscores using different methods based on original CTT scores.}
\usage{
CTTsub(test.data, method = "Haberman")
}
\arguments{
\item{test.data}{A list that contains item responses of all subtests and the entire test, 
which can be obtained using function 'data.prep'.}

\item{method}{Subscore estimation methods. method="Haberman" (by default) represents 
the three methods proposed by Haberman (2008) <doi:10.3102/1076998607302636>. 
method="Wainer" represents Wainer's augmented method.}
}
\value{
\item{summary}{Summary of estimated subscores (e.g., mean, sd).}
\item{PRMSE}{(a) PRMSE values of estimated subscores (for Haberman's methods only).(b) Decisions on whether subscores have added 
value - added.value.s (or added.value.sx) = 1 means subscore.s (or subscore.sx) has added value, and  added.value.s (or added.value.sx) = 0
vice versa.}
\item{PRMSE.test}{All information in PRMSE plus results of hypopthesis testing based on Sinharay (2019) <doi:10.3102/1076998618788862>.}
\item{subscore.original}{Original subscores and total score.}
\item{estimated.subscores}{Subscores computed using selected method. Three sets of subscores will be returned if method = "Haberman".}
}
\description{
This function estimates true subscores using methods introduced in studies of
Haberman (2008) <doi:10.3102/1076998607302636> and Wainer et al. (2001) <doi:10.4324/9781410604729>. 
Hypothesis tests (i.e., Olkin' Z,Williams's t, and Hedges-Olkin's Z) are used to determine 
whether a subscore or an augmented subscore has added value. Codes for the hypothesis tests are from Sinharay (2019) 
<doi: 10.3102/1076998618788862>.
}
\examples{
# Transferring original scored data to a list format
# that can be used in other functions.
test.data<-data.prep(scored.data,c(3,15,15,20),
                     c("Algebra","Geometry","Measurement", "Math"))
#----------------------------------------------
# Estimating subscores using Haberman's methods       
CTTsub(test.data,method="Haberman") # Estimating subscores using Haberman's methods 

# Obtaining original correlation for the three methods
CTTsub(test.data,method="Haberman")$Correlation  

# Obtaining disattenuated correlation for the three methods
CTTsub(test.data,method="Haberman")$Disattenuated.correlation  

# Obtaining PRMSEs for the three methods
CTTsub(test.data,method="Haberman")$PRMSE  

# Obtaining descriptive statistics summary for estimated subscores  
CTTsub(test.data,method="Haberman")$summary 

# Obtaining raw subscores  
CTTsub(test.data,method="Haberman")$subscore.original  

# Obtaining subscores that are estimated as a function of the observed subscores 
CTTsub(test.data,method="Haberman")$subscore.s 

# Obtaining subscores that are estimated as a function of the observed total score 
CTTsub(test.data,method="Haberman")$subscore.x  

# Obtaining subscores that are estimated as a function of 
# both the observed subscores and the observed total score.
CTTsub(test.data,method="Haberman")$subscore.sx  

#-------------------------------------------      
# Estimating subscores using Wainer's method
CTTsub(test.data,method="Wainer") 
       
# Obtaining descriptive statistics summary for subscores
CTTsub(test.data,method="Wainer")$summary   

# Obtaining original subscores
CTTsub(test.data,method="Wainer")$subscore.original 

# Obtaining subscores that are estimated using Wainer's augmentation method  
CTTsub(test.data,method="Wainer")$subscore.augmented  
        
}
\references{
{
Haberman, S. J. (2008). 
"When can subscores have value?."
 Journal of Educational and Behavioral Statistics, 33(2), 204-229. doi:10.3102/1076998607302636. 
}

{
Sinharay, S. (2019). 
"Added Value of Subscores and Hypothesis Testing."
 Journal of Educational and Behavioral Statistics, 44(1), 25-44. doi:10.3102/1076998618788862.
}

{
Wainer, H., Vevea, J., Camacho, F., Reeve, R., Rosa, K., Nelson, L., Swygert, K., & Thissen, D. (2001). 
"Augmented scores - "Borrowing strength" to compute scores based on small numbers of items."
 In Thissen, D. & Wainer, H. (Eds.), Test scoring (pp.343 - 387). Mahwah, NJ: Lawrence Erlbaum Associates, Inc. 
 doi:10.4324/9781410604729.
}
}
