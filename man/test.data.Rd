\name{test.data}
\alias{test.data}
\docType{data}
\title{
A list of objects that include both test information and subscores.
}
\description{
This list consists of four objects. The first three objects are item responses of the three subtests. 
The fourth object is the scored data of the whoel test.
}
\usage{data("test.data")}
\format{
  The format is:

List of 4

 $ subtest.1 :'data.frame':	150 obs. of  15 variables:

 $ subtest.2 :'data.frame':	150 obs. of  15 variables:

 $ subtest.3 :'data.frame':	150 obs. of  20 variables:

 $ total.test:'data.frame':	150 obs. of  50 variables:
  
}
\details{
Object 1: 150 participants and 15 items; Object 2: 150 participants and 15 items.
Object 3: 150 participants and 20 items; Object 4: 150 participants and 20 items.
}
\examples{
data(test.data)
## maybe str(test.data) ; plot(test.data) ...
}
\keyword{datasets}
