library(randomForest)
library(forestFloor)
library(rfFC)
library(trimTrees)
obs=2000
vars=7
X = data.frame(1:obs)
for(i in 1:vars) X[,i] = sample(1:7,obs,replace=T)
Y = apply(X,1,function(x) mean(x))
Y[Y>mean(Y)]=42
Y[Y<mean(Y)]=666
hist(Y)
Y=as.factor(Y)

for(i in 1:vars) X[,i]=as.factor(X[,i])
rfo = cinbag(X,Y,keep.inbag=T,replace=T,ntree=800)
ff = forestFloor2(rfo,X,T)
plot(ff)
plot(apply(ff$FCmatrix,1,sum),Y)
cor(rfo$predicted,apply(ff$FCmatrix,1,sum))

featureContributions
li = getLocalIncrements(rfo,X)
fc = featureContributions(rfo,li,X)
fcpred = apply(fc,1,sum) 
plot(predict(rfo,X),fcpred,cex=0.1)
fcpred

