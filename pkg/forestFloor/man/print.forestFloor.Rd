\name{print.forestFloor}
\alias{print.forestFloor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
print summary of forestFloor.Object
}
\description{
summary comprise N_fold I_repeated cross validation statistics and aggragated importance of variables.
a forestFloor_Object is of class forestFloor and is the output of the forestFloor function.
}
\usage{
print(forestfloor_Object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{forestFloor_Object}{
%%     ~~Describe \code{RFCV_object} here~~
}
}
\details{
prints short help text for usage of a forestFloor_object
}

\author{
Soren Havelund Welling
}


\examples{
#simulate data
obs=1000 \cr
vars = 6 \cr
X = data.frame(replicate(vars,rnorm(obs)))\cr 
Y = with(X, X1^2 + sin(X2*pi) + 2 * X3 * X4 + 0.5 * rnorm(obs))\cr 

#grow a forest, remeber to include inbag
library(randomForest)
rfo=randomForest(X,Y,keep.inbag=TRUE)\cr

#compute topology
ff = forestFloor(rfo,X)\cr

#print forestFloor
print(ff) \cr


}