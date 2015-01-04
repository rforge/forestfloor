\name{forestFloor-package}
\docType{package}
\title{
forestFloor: visualize the randomForest topology
}
\description{
forrestFloor visualizes cross-validated topology-maps of randomForests(RF).
  Package enables users to understand a non-linear, regression problem or a binary
  classification problem through RF. In all, this package is intended to provide a fast
  overview of dynamics within a given system of interest, allowing the user to decide
  for apropiate further modeling  maybe within a classical statistical framework or to stay
  within the RF-modeling and look deep into  the alluring topology of correlations and
  local interactions.
}
\details{
\tabular{ll}{
Package: \tab forestFloor\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-12-07\cr
License: \tab GPL-2\cr
}
}

\author{
Søren Havelund Welling
}
\references{
Interpretation of QSAR Models Based on Random Forest Methods, http://dx.doi.org/10.1002/minf.201000173 \cr
Interpreting random forest classification models using a feature contribution method, http://arxiv.org/abs/1312.1121 \cr
}

\examples{
#simulate data
obs=1000 
vars = 6 
X = data.frame(replicate(vars,rnorm(obs))) 
Y = with(X, X1^2 + sin(X2*pi) + 2 * X3 * X4 + 0.5 * rnorm(obs)) 

#grow a forest, remeber to include inbag
rfo=randomForest(X,Y,keep.inbag=TRUE)

#compute topology
ff = forestFloor2(rfo,X)

#print forestFloor
print(ff) 

#plot partial functions of most important variables first
plot(ff,plot_seq=NULL,colour_by="dummy, no defined color-template",col_axis=1) 

#Non interacting functions are well displayed, whereas X3 and X4 are not
#by applying different colourgradient, interactions reveal themself 
plot(ff,plot_seq=NULL,colour_by=3,col_axis=1) 

#in 3D the interaction between X3 and X reveals itself completely
show3d(ff,x_cols=3,y_cols=4,z_cols=3:4) 

#although no interaction, a joined additive effect of X1 and X2
#can also be informative to display in 3D
plot(ff,plot_seq=NULL,colour_by=4,col_axis=2) #use plot first to define colours
show3d(ff,x_cols=1,y_cols=2,z_cols=1:2)
}