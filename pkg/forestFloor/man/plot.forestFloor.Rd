\name{plot.forestFloor}
\alias{plot.forestFloor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Metod: plot.forestFloor
}
\description{
Method to plot an object of forestFloor-class. Plot partial feature contributions of the most important variables. Colour gradients can be applied two show possible interactions.
}
\usage{
\method{plot}{forestFloor}(x, colour_by = 1, col_axis = 1, plot_seq = NULL, alpha = "auto", limitY = TRUE, order_by_importance = T, external.col = NULL, 
    cropXaxes = NULL, crop_limit = 4, ...) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
foretFloor-object, also abbrivated ff..
Computed topology of randomForest-model, the output from the forestFloor function \cr
includes also X and Y and importance data
}

  \item{colour_by}{
different colour method can be applied: \cr
input is a integer or "PCA" or "top2" \cr
a number will point to one specific variable, a gradient from red to green to blue is applied. colour_by = 2 will refer to the second of multiple plots. Which plot is the second is influenced by 'plot_seq' and 'order_by_importance' see below. 
"PCA" will apply importance weighted PCA, to reduce the system to 3 linear principal components which is converted into RGB.color. Similarly coloured points are likely to be linearly close.
"top2" applies a two.way colouring gradient by the most and second most important variable
}

\item{col_axis}{
a number \cr
=1 gradient is computed from feature space or X-axis \cr
=2 gradient is compputed from feature contribution space y-axis \cr

... if colour_by = 2 & col_axis = 2 a vertical gradient will be applied to second plot (and to others correspondingly). \cr
... if colour_by = "PCA" & col_axis = 1, three feature space PCA components will colour the observations \cr
... if colour_by = "PCA" & col_axis = 2, three feature contributions PCA components will colour the observations \cr
}

\item{plot_seq}{
a numeric vector describing which variables and in what sequence to plot,
ordered by importance as default, order_by_importance = F then by feature/coloumn order of training data.
}

\item{alpha}{
a number within [0;1] or the method "auto" 
transparency of points,
"auto" the more points the more transparrent to avoid overplotting
otherwise manually set as a number [0;1] where 0 is completely transparrent
}

\item{limitY}{
TRUE/FLASE, constrain all Yaxis to same limits to ensure relevance of low importance features is not overinterpreted
}

\item{order_by_importance}{
TRUE / FALSE
should plotting and plot_seq be ordered after importance. Most important feature plot first(TRUE)
}

\item{external.col}{
a colour vector being sent to col input of genric plot function plot.
This port allows a user to define new colour gradients. e.g.
ext.col.vec = rgb(red=box.outliers(Y),green=0.2,blue=0.2,alpha=0.4) \cr
Where Y is the numeric training target, box.outliers() remove extreme observations and normalize, rgb(). ext.col.vec should in any case have same length as Y.
}

\item{cropXaxes}{
a vector of indice numbers of which zooming of x.axis should look away from outliers
}

\item{crop_limit}{
a number often between 1.5 and 5, referring limit in std.devs from the mean defining outliers
if limit = 2, above selected plots will zoom to +/- 2 std.dev of the respective features.
}

\item{...}{
... other arguments passed to generic plot functions
}

}
\details{
The method plot.forestFloor visualizes partial plots of the most important variables first. Partial plotting is also available in randomForest package. But such plots are single thin lines and do not answer the questions;
Is this partial function(PF) a fair generalization or subject to global or local interactions, which cannot be captured in one_way partial plots? \cr
Do sufficently datapoints describe all areas of the PF, or could some areas of the partial plot have a small observation density? \cr
Each partial plot should be viewed as perspectives of the multidimensional toplogy. With multiple perspectives with matching colour gradients and occasional 3D-plots, a fairly complex high dimensional topology can be understood.\cr
Start inspecting the most important but diffuse partial function and recolour observations by this variables physical value. If systematic colourings appear in other partial plots, these two variables/feature contributions are related. If the colour-gradiant most appear horisontal, the variables are correlated. If the colour pattarns are stacking vertical, these feature contributions is interacting in this region. This mean the random forest topology applies an interaction between these two variables. Taking such 2Way interaction further to 3D can comfirm such interaction(see show3D). Looking back to 2D coloured plots, the user can now recognize the shadow of 3D structure. Infact the colouring reminds of that structure. There is no easy way to plot even higher level interactions. But maybe the user can use 3D and colours to percieve 4D and 5D aswell, as we here can use 2D+colours to depict some higher interactions. \cr

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

\keyword{ multivariate }
\keyword{ models }
\keyword{ nonlinear }
\keyword{ robust }
