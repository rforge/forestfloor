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
show3d(x, order_by_importance = FALSE, which_matrices = c("X", 
    "X", "FCmatrix"), x_cols = 1, y_cols = 2, z_cols = c(1:2), 
    plot.surface = TRUE, grid.lines = 30, k = 5, alpha.surf = 0.4, 
    alpha.obs = 0.4, size.obs = 3, z_scale = 0.7, knnBag = 20, 
    bag.ratio = 0.5, avoidFreeType = TRUE, ...) 
%- maybe also 'usage' for other objects documented here.
}
\arguments{
  
  \item{x}{
  also abbrivated ff. 
  Computed topology of randomForest-model, the output object from the forestFloor function of class forestFloor
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
vector of one or multiple col indices for z_axis. By default col indices refers training feature contribution matrix ordered by importance. z_cols would most times hold both or some of the indice numbers from x_cols and y_cols.
  }
  \item{plot.surface}{
TRUE/FALSE
  }
  \item{grid.lines}{
integer 2 to many. How many grid.lines in surface response. 30 grid.lines = 30^2 = 900 anchorpoints
  }
  \item{k}{
integer from 1 to many. k amount of nearest neighbors to estimate response of anchor point in surface. Higher values will smooth surface and make it biased. Lower values will be increase noise.
  }
  \item{alpha.surf}{
a number within [0;1]
transparency of surface, lower values more transparent
  }
  \item{alpha.obs}{
a number within [0;1]
transparency of observations, lower values more transparent
 }
  \item{size.obs}{
A postive number as the size of observations.
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
  
  \item{avoidFreeType}{
  avoidFreeType pass specificly to rgl function not to use vector.based labels. \cr
  freeType make labels look better but must installed externally in operating system and it could take some time.
  }
  
    \item{...}{
  ... other arguments passed to generic rgl-functions rgl::pers3d and rgl::points3d  or rgl::plot3d when plot_surface=FALSE.
  }
  
}
\details{
show3d allows to plot data from ff, the forestFloor-object.
Most users will only need to change following:
x_cols, y_cols the featurespace axis \cr
z_cols the feature contributions axis \cr
...and sometimes these arguments: \cr
k the smoothness of the surface\cr
surface was very simple to implement, in long run would true gaussian weighting paremeter maybe be advantageous.
}

\value{
no value

}

\author{
Soren Havelund Welling
}

\examples{
#simulate data
obs=1000 
vars = 6 
X = data.frame(replicate(vars,rnorm(obs))) 
Y = with(X, X1^2 + sin(X2*pi) + 2 * X3 * X4 + 0.5 * rnorm(obs)) 

#grow a forest, remeber to include inbag
rfo=randomForest::randomForest(X,Y,keep.inbag=TRUE)

#compute topology
ff = forestFloor(rfo,X)

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
