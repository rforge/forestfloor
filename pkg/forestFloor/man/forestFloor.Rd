\name{forestFloor}
\alias{forestFloor_f}
\alias{forestFloor}

\title{
ForestFloor: Visualize topologies of randomForest model
}
\description{
The function forestFloor computes a feature contribution matrix from a randomForest model-fit and outputs a forestFloor S3 class object, including importance and the orginal training set. The output object is the basis all visualizations.

}
\usage{
ff = forestFloor(rfo,X,calc_np=FALSE)
}

\arguments{
\item{rfo}{
rfo, random forest object is the output from randomForest::randomForest or cinbag::trimTrees \cr
for regression use: rfo = randomForest(X,Y,keep.inbag=T,importance=T) \cr
for binary classification use: rfo = cinbag(X,Y,keep.inbag=T,keep.forest=T,importance=T) \cr
Formula interface is not supported. Y is a factor of two levels or numeric vector or n_row elements. Multi-class clasification is not supported, try rfFC-package instead see references.
}

\item{X}{
data.frame of input variables, numeric(continnous), descrete(treated as continous) or factors(categoric).
n_rows observations and n_columns features
X MUST be the same data.frame as used to train the random forest, see above item.
}

\item{calc_np}{
  calculate Node Predictions(TRUE) or reuse information from rfo(FALSE)?

  slightly faster when FALSE for regression \cr
  MUST be TRUE for binary classification \cr

  Node predictions, the average target value of inbag samples in any terminal or intermediary node of a random forest are already calculated for regression and are placed in rfo$forest$nodepred. Node predictions for binaray classification are the fraction of class 1 (out of 2) in any node of a random forest and are not calculated in advance.
}
}

\details{
  forestFloor computes feature contributions for random forest regression as suggest by Kuz'min et al, and for binaray classification as suggested by Palczewska et al. Feature contributions is the sums over all local increments for each observation for each feature divided by the number of trees. A local increment is the change of node prediction for given observation in one node being split to a subnode by a given feature. forestFloor use inbag samples to calculate local increments, but only sum local increments over out-of-bag samples divided with OOBtimes. OOBtimes is the number of times a given observation have out-of-bag which normally is ~ trees / 3. This implementation, can be said to yield cross-validated feature contributions. In practices this lowers the leaverage of any observation to the feature contributions of this observation. Hereby becomes the visulization less noisy. In systems with low or no noise, this implementation have no particular advantage.
}

\value{
  the forestFloor function outputs an object of class "forestFloor" with following elements:

\item{X}{
  a copy of the training data or feature space matrix/data.frame, X. The copy is passed from the input of this function. X is used in all visualization to expand the feature contributions over the features of which they were recorded.
}

\item{Y}{
  a copy of the target vector, Y.
}

\item{importance}{
  The gini-importance or permutation-importance a.k.a varaiable importance of the random forest object \cr
  if rfo=randomForest(X,Y,importance=FALSE), gini-importance is used. \cr
  gini-importance is less reproducible and more biased. The extra time used to compute permutation importance is     negliable. 
}

\item{imp_ind}{
  imp_ind, the importance indices is the order to sort the features by descending importance. imp_ind is used by plotting functions to present must relevant feature contributions first. If using gini-importance, the order of plots is more random and will favor continous variables. The plots themselves will not differ.
}

\item{FC_matrix}{
  feature contributions in a matrix. \cr
  n_row observations and n_column features - same dimensions as X.
}

}

\references{
Interpretation of QSAR Models Based on Random Forest Methods, http://dx.doi.org/10.1002/minf.201000173 \cr
Interpreting random forest classification models using a feature contribution method, http://arxiv.org/abs/1312.1121 \cr
}

\author{
Søren Havelund Welling, November 2014
}

\note{
this version 
}

\seealso{
\code{\link{plot.forestFloor}}, \code{\link{show3d}}
}
\examples{

#simulate data
obs=1000 \cr
vars = 6 \cr
X = data.frame(replicate(vars,rnorm(obs)))\cr 
Y = with(X, X1^2 + sin(X2*pi) + 2 * X3 * X4 + 0.5 * rnorm(obs))\cr 

#grow a forest, remeber to include inbag
rfo=randomForest(X,Y,keep.inbag = TRUE)\cr

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

\keyword{ multivariate }
\keyword{ models }
\keyword{ nonlinear }
\keyword{ robust }