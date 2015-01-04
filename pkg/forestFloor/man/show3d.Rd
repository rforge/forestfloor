\name{show3d}
\alias{show3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
make 3d-plot of forestFloor topology
}
\description{
Each axis is provided one variable or feature contribution, a z-axis surface response layer will be estimated and plotted alongside the datapoints. Latests 2D colour gradient is fetched from object global.col in .GlobalEnv. show3d() uses openGL rendering through the rgl-package. Image will emerge in external graphical device and can controlled by mouse pointer.

}
\usage{
show3d(RFCV_object, order_by_importance = F, which_matrices = c("X", "X", "FCsCV"), x_cols = 2, y_cols = 11, z_cols = c(2), plot.surface = T, grid.lines = 30, k = 5, alpha = 0.4, z_scale = 0.7, knnBag = 20, bag.ratio = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{forestFloor_object}{
  Computed topology of randomForest-model, the output from the forestFloor function of class forestFloor
%%     ~~Describe \code{RFCV_object} here~~
  }
  \item{order_by_importance}{
   = FALSE: indices refer to variable order in original training data. (good for scripting) \cr
   = TRUE:  indices refer to importanced ordered. (good for user interaction)
  }
  \item{which_matrices}{
  a char vector of length 3, with names of object data should be fetched. As default is 
  c("X", "X", "FCsCV"); indices of X_cols and y_cols will be fetched from training set matrix and z_cols from feature contribution matrix. Any matrix in forstFloor.object can theoretically accessesible for plotting. This setting is rarely changed.
  
  }
  \item{x_cols}{
vector of one or multiple col indices for x_axis. By default col indices refers training variable matrix ordered by importance.

  }
  \item{y_cols}{
vector of one or multiple col indices for x_axis. By default col indices refers training variable matrix ordered by importance.

  }
  \item{z_cols}{
vector of one or multiple col indices for z_axis. By default col indices refers training feature contribution matrix ordered by importance.
  }
  \item{plot.surface}{
TRUE/FALSE
  }
  \item{grid.lines}{
integer 2 to many. How many grid.lines in surface response. 30 grid.lines = 900 anchorpoints
  }
  \item{k}{
interger from 1 to many. k amount of nearest neighbors to average response of anchor point. Higher values will smooth surface.
  }
  \item{alpha}{
transparency of surface, lower values more transparent
  }
  \item{z_scale}{
numeric, near zero to larger number. Steapness of surface. Between 0.1 and 10 are reasonable settings. Paremeter is passed directly to rgl-functions.
  }
  \item{knnBag}{
integer, 1 to many. How many times should kNN be bootstrapped and estimated for each anchorpoint. High values increase reproducibility of surface rarely more than 20 is needed.
  }
  \item{bag.ratio}{
numeric number between near zero and one. The ratio of used observations over all observations used in each bootstrapped surface estimation. Values less than one, will increase amount above k nearest neighbors used in estimation of each anchorpoint. Nearest neighbors will participate in most estimations and have a high weighting, more distant neighbors will participate gradually more rarely and yield a lower weighting. bag.ratio is therefore a non-metric pseudo gaussian weighting paremeter.
  }
}
\details{
was very simple to implement, in long run would true gaussian weighting paremeter maybe be advantageous. Computing times of bagge nearest neighbor is theoretically poor but practically instancetaously.
}

\value{
no value
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Søren Havelund Welling
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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

#plot partial functions of most important variables first
plot(ff,plot_seq=NULL,colour_by="dummy, no defined color-template",col_axis=1) \cr

#Non interacting functions are well displayed, whereas X3 and X4 are not
#by applying different colourgradient, interactions reveal themself 
plot(ff,plot_seq=NULL,colour_by=3,col_axis=1) \cr

#in 3D the interaction between X3 and X reveals itself completely
show3d(ff,x_cols=3,y_cols=4,z_cols=3:4) \cr

#although no interaction, a joined additive effect of X1 and X2
#can also be informative to display in 3D
plot(ff,plot_seq=NULL,colour_by=4,col_axis=2) #use plot first to define colours \cr
show3d(ff,x_cols=1,y_cols=2,z_cols=1:2) \cr

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
