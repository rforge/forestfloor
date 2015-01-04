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
\method{plot}{forestFloor}(ff, colour_by = 1, col_axis = 1, plot_seq = NULL, alpha = "auto", 
    limitX = FALSE, limitY = TRUE, order_by_importance = T, external.col = NULL, 
    cropXaxes = NULL, crop_limit = 4, ...) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ff}{
foretFloor-object
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
}
\details{
The method plot.forestFloor visualizes partial plots of the most important variables first. Partial plotting is also available in randomForest package. But such plots are single thin lines and do not answer the questions;
Is this partial function(PF) a fair generalization or subject to global or local interactions, which cannot be captured in one_way partial plots? \cr
Do sufficently datapoints describe all areas of the PF, or could some areas of the partial plot have a small observation density? \cr
Each partial plot should be viewed as perspectives of the multidimensional toplogy. With multiple perspectives with matching colour gradients and occasional 3D-plots, a fairly complex high dimensional topology can be understood.\cr
Start inspecting the most important but diffuse partial function and recolour observations by this variables physical value. If systematic colourings appear in other partial plots, these two variables/feature contributions are related. If the colour-gradiant most appear horisontal, the variables are correlated. If the colour pattarns are stacking vertical, these feature contributions is interacting in this region. This mean the random forest topology applies an interaction between these two variables. Taking such 2Way interaction further to 3D can comfirm such interaction(see show3D). Looking back to 2D coloured plots, the user can now recognize the shadow of 3D structure. Infact the colouring reminds of that structure. There is no easy way to plot even higher level interactions. But maybe the user can use 3D and colours to percieve 4D and 5D aswell, as we here can use 2D+colours to depict some higher interactions.

Why not just plot the raw data as is, in various ways? Because even when noise component is low, the plotted 2D shadow of multiple non-linear interactions can look merely random. randomForest provides interpretation of what there is the signal and of individual contributions/components related to each variable. What would otherwise become a messy indistinguishable cloud of interactions and additive effects and noise.\cr
system = signal + noise = additive component + interacting components + noise \cr
This pseudo-equation explains how a high dimensional system can be segrated into smaller parts which in turn can be visulized. The random forest provides an interpretation of what is signal and feature contributions of the individual feature related components. Adittive components can be plotted as is by this function. Interacting components must be plotted in the context of two or more features, see \code{\link{show3d}}

The root of computation of the partial interpretations, is that treesbased models, recursively process the data univariately - one dicision by one variable at the time. This leads to a multitude of interpretations where randomForest is the average of many trees and feature contributions is the sampled topology of this fit to the data. Secondly random forest is robust and thereby removing a good deal of the noise component. Random trends will across many trees average eachother out, wheres reproducible trends will stay. Furthermore default paremeter-set of randomForest-algorithm is stable and do not need much optimising. Other methods may later supercede in predictions performance. If an underlying system is near linear, GLM will always perform better. But before blindly choosing between e.g. RF, SVM and GLM solely by prediction performance, an idea of a expected system topology should be highly appreciated. For an archtypical 'conservative classical statistician', this package would be funpark trip, only inspirring to what models should tested out to be significant. Could any of the variables be transformed to become linaer? What interactions-terms should be included?

Interpreting partial feature/variable contributions can be an inspirring or alleuring tool. Plausible causality links comes to mind when interactions are mapped. Remember to also concider reverse and external causality or that sampling was not truely indpendent, e.g. Simpson Paradox. Furthermore oberserved interactions can be due to general statistical issues, such as high colinearity, low sample size and non-independant sampling and more. Often some usaul supects can be ruled out as implausible, while others could be similar plausible and only further testing could possibly tell. Remember these topology maps are created by just another oppurtunistic algorithms, only trying to please its loss function. When data is sparse, randomForest tends to yield topologies with few interactions and with steap sigmoid partial functions. This does not mean the underlying system is exactly as such, just that no more complex interpretations was stable and reproducible across the many treees. Furthermore tree models do not extrapolate outside featurespace of trainingset. Any datapoint outside will be predicted as the most resembling datapoints inside the feature space. The amplitude of suggested partialfunctions is lowered in the near boarder of feature space until et becomes a flatline. This make linear effects look more or less sigmoid. This effect is very well displayed with simulated sinfunctions which have large amplitudes in the center of feature space and small gradually dimishing to no amplitude close to the border of feature space. Also normal distribted variables increase this effect of soft border in the feature space. In high dimensions featurespace of uncorrealted variales will have a low density, and if the variables also are normal distributed, the density of datapoints at the borders of features space will be even more sparse. Personally I think that any data-driven model should return to conservative estimates in low density areas of feature space, which is what randomForest in generally do. Nonetheless, as randomForest is not completely robust, sometimes some complex topologies can emerge by chance in a low density areas of feature space. Here it is important to look at the amount data point driving this specifc topology. With a given prior belief is this a plausible and/or probable topology? On th otherhand when firm topologies are found and in areas of many data points, the rule of large numbers must apply, and this should emphasised more closely.
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
Søren Havelund Welling, 2014-15
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
obs=1000
vars = 6 
X = data.frame(replicate(vars,rnorm(obs)))
Y = with(X, X1^2 + sin(X2*pi) + 2 * X3 * X4 + 0.5 * rnorm(obs))

#grow a forest, remeber to include inbag
rfo=randomForest(X,Y,keep.inbag=TRUE)

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
