if(!interactive()) Sys.setenv(RGL_USE_NULL=TRUE) #disable RGL for headless machines
library(forestFloor)
library(randomForest)
#simulate data
obs=2000
vars = 6 

X = data.frame(replicate(vars,rnorm(obs)))
Y = with(X, X1^2 + sin(X2*pi) + 2 * X3 * X4 + .5 * rnorm(obs))


#grow a forest, remeber to include inbag
rfTest42=randomForest(X,Y,keep.inbag = TRUE,sampsize=1000,ntree=500)

#compute topology
ffTest42 = forestFloor(rfTest42,X)


#print forestFloor
print(ffTest42) 

#plot partial functions of most important variables first
plot(ffTest42,orderByImportance=TRUE) 

#Non interacting functions are well displayed, whereas X3 and X4 are not
#by applying different colourgradient, interactions reveal themself 
#also a k-nearest neighbor fit is applied to evaluate goodness of fit
Col=fcol(ffTest42,3,orderByImportance=FALSE)
plot(ffTest42,col=Col,plot_GOF=TRUE,speed=T) 

#if ever needed, k-nearest neighbor parameters for goodness-of-fit can be access through convolute_ff
#a new fit will be calculated and added to forstFloor object as ffTest42$FCfit
ffTest42 = convolute_ff(ffTest42,userArgs.kknn=alist(kernel="epanechnikov",kmax=5))
plot(ffTest42,col=Col,plot_GOF=TRUE)

#in 3D the interaction between X3 and X reveals itself completely
show3d(ffTest42,3:4,col=Col,plot.rgl=list(size=5),orderByImportance=FALSE) 

#although no interaction, a joined additive effect of X1 and X2
#colour by FC-component FC1 and FC2 summed
Col = fcol(ffTest42,1:2,orderByImportance=FALSE,X.m=FALSE,RGB=TRUE)
plot(ffTest42,col=Col) 
show3d(ffTest42,1:2,col=Col,plot.rgl=list(size=5),orderByImportance=FALSE) 

#...or two-way gradient is formed from FC-component X1 and X2.
Col = fcol(ffTest42,1:2,orderByImportance=FALSE,X.matrix=TRUE,alpha=0.8) 
plot(ffTest42,col=Col) 
show3d(ffTest42,1:2,col=Col,plot.rgl=list(size=5),orderByImportance=FALSE,plot_GOF=T)

