
# Methods:
#m1 print output
print.forestFloor = function(x,...) {
 cat("this is a forestFloor('x') object \n
this object can be plotted in 2D with plot(x), see help(plot.forestFloor) \n
this object can be plotted in 3D with show3d(x), see help(show3d) \n
\n
x contains following internal elements: \n ",with(x,ls()))
}

#m2 plot output
plot.forestFloor = function(x,
                            colour_by=1,
                            col_axis = 1,
                            plot_seq=NULL,
                            alpha="auto",
                            limitY=TRUE,
                            order_by_importance=T,
                            external.col=NULL,
                            cropXaxes=NULL,
                            crop_limit=4,
                            ...)
  {
  
  pars = par(no.readonly = TRUE) #save previous graphical par(emeters)
  par(mar=c(2,2,1,1),cex=.5) #changing par, narrowing plot margins, smaller points
  
  #short for phys.val and feature contribution in object
  X = x$X
  FCs = x$FCmatrix
  
  #Auto setting transparancy variable. The more obs, the more transparrency
  if(alpha=="auto") alpha = min(max(400/dim(X)[1],0.2),1)
  
  #If now sequnce, choosing to plot first 18 variables
  if(is.null(plot_seq)) plot_seq = 1:min(dim(X)[2],24)
  
  #make catogorical features numeric, save jitter.template
  jitter.template=rep(FALSE,dim(X)[2]) #list of what features are catagorical
  as.numeric.factor <- function(x) {match(x,levels(x))}
  for(i in 1:dim(X)[2]) {
    if(is.factor(X[,i])) {
      jitter.template[i]=T
      this.fac=as.numeric.factor(X[,i])
      X[,i] = this.fac
    }
    if(is.character(X[,i])) X[,i] = as.numeric(X[,i])
  } 
  
  ##get dimensions of plots
  n.plots = min(dim(X)[2],length(plot_seq))
  plotdims.y = min(ceiling(n.plots/3),5)
  plotdims.x = min(3 , n.plots)
  par(mfrow=c(plotdims.y,plotdims.x))
  
  ##get importance for plotting
  imp = x$importance     #fetch importance
  imp.ind = x$imp_ind    #fetch importance ranking/indices
  
  
  #set default.colour
  colours = "black"
  if(!is.null(external.col)) {
    colours = external.col
  } else {
    
  #following is multiple colour rules:
  
  if(colour_by=="PCA") {
  #make importance scaled PCA of X
    if(col_axis==1) {
      pca.X  = prcomp(scale(X)*t(replicate(dim(X)[1],imp)))
    } else {
      pca.X  = prcomp(scale(FCs)*t(replicate(dim(FCs)[1],imp)))
    }
      PC123 = pca.X$x[,1:3] #fetch 3 first principal components

  #apply box.outliers to PCA components, to avoid extreme colour leverage of outliers
  PC123.box = apply(PC123,2,box.outliers)
  #change into colours, with transparancy
  colours   = apply(PC123.box,1,function(x) rgb(x[1],x[2],x[3],alpha=alpha))
  }
  
  #colour by top2 variables
  if(colour_by=="top2") {
    if(col_axis==1) {
      sX = apply(X[,imp.ind[1:3]],2,box.outliers)
    }else{
      sX = apply(FCs[,imp.ind[1:3]],2,box.outliers) 
    }
    nX = apply(sX,2, function(x) (x-min(x))/(max(x)-min(x)))
    nX = apply(nX,2, function(x) (x+.5) / 1.5)
    colours = apply(nX,1,function(x) rgb(x[1],x[2],.1,alpha=alpha))
  }
  
  ## colour by specific variable
  if(is.numeric(colour_by)) {
    if(col_axis==1) {
      sX = apply(X[,imp.ind[rep(colour_by,3)]],2,box.outliers)
    } else {
      sX = apply(FCs[,imp.ind[rep(colour_by,3)]],2,box.outliers)
    }
    
    nX = apply(sX,2, function(x) (x-min(x))/(max(x)-min(x)))
    colours = apply(nX,1,function(x) rgb(x[1]^3,1-x[1]^3-(1-x[3])^3,(1-x[3])^3,alpha=alpha))
    #rainbow colouring
    #nx = X[,imp.ind[colour_by]]
    #colours = rainbow(length(nx),start=0.2,alpha=alpha)[match(1:length(nx),sort(nx,index.return=T)$ix)]
  }
  
  }
  
  #Save this colouring globally, for later 3D plotting
  if(exists("forestFloor_graphics.env",envir=.GlobalEnv)) {
    assign("obs.indv.colours",colours,envir=forestFloor_graphics.env)
  } else {
    local({forestFloor_graphics.env <- new.env()},envir=.GlobalEnv)
    assign("obs.indv.colours",colours,envir=forestFloor_graphics.env)
    #goto global env, make graphics.env, place colours here
  }
  
  
  
  
  ##plot the n.plots most important variables
  Xsd = 0:1 #initialize Xsd
  if(!order_by_importance) imp.ind=sort(imp.ind) #optinal removal of importance order
  for(i in plot_seq) {
    
    if(i %in% cropXaxes && !is.null(cropXaxes)) {
      limitX = T
      Xsd = box.outliers(as.numeric(X[,imp.ind[i]],2),limit=crop_limit,normalize=F)
    } else {
      limitX = F
    }
    
    plot(
      data.frame( # data to plot
        physical.value        = jitter(X[,imp.ind[i]],factor=jitter.template[imp.ind[i]]*2),
        partial.contribution  = FCs[,imp.ind[i]]
      ),
      main = names(imp)[imp.ind[i]],
      col = colours,  #colours are fetched from forestFloor_graphics.env
      ylim = list(NULL,range(FCs))[[limitY+1]], #same Yaxis if limitY == TRUE
      xlim = list(NULL,range(Xsd))[[limitX+1]],
      ...
    )
  }
  #print(pars$pin)
  par(pars)
}

#f1 3d show function
show3d = function(x,
                  order_by_importance=FALSE,
                  which_matrices=c("X","X","FCmatrix"),
                  x_cols=1,y_cols=2,z_cols=c(1:2),
                  plot.surface=TRUE,
                  grid.lines=30,
                  k=5,
                  alpha.surf=.4,
                  alpha.obs=.4,
                  size.obs=3,
                  z_scale=.7,
                  knnBag=20,
                  bag.ratio=0.5,
                  avoidFreeType = T,
                  ...) {
  
  #retrieve labels for plotting
  xyzlab = with(x,{
    lab=c()
    if(!order_by_importance) imp_ind = sort(imp_ind)
    lab$x =       names(get(which_matrices[1]))[imp_ind[x_cols]]
    lab$y =       names(get(which_matrices[2]))[imp_ind[y_cols]]
    lab$z =       names(get(which_matrices[3]))[imp_ind[z_cols]]
    return(lab)
  })
  
  
  #apply funciton later requires two inputs pr axis, if one that is duplicated into two
  cols.fix = function(i) if(length(i)<2){i=c(i,i)}else{i}
  x_cols = cols.fix(x_cols)
  y_cols = cols.fix(y_cols)
  z_cols = cols.fix(z_cols)
  
  
  axisval = with(x,{
    axisval=c()
    if(order_by_importance) {
      imp.ind = imp_ind
      ind_by_imp = function(cols,imp.ind) imp.ind[cols]
      x_cols = ind_by_imp(x_cols,imp.ind)
      y_cols = ind_by_imp(y_cols,imp.ind)
      z_cols = ind_by_imp(z_cols,imp.ind)
    }
    
    if(!which_matrices[1]=="FCmatrix") {
      axisval$x =          get(which_matrices[1])[,x_cols[1]]
    }else{
      axisval$x =    apply(get(which_matrices[1])[,x_cols],1,mean)
    }
    if(!which_matrices[2]=="FCmatrix") {
      axisval$y       =    get(which_matrices[2])[,y_cols[1]]
    }else{
      axisval$y =    apply(get(which_matrices[2])[,y_cols],1,mean)
    }
    if(!which_matrices[3]=="FCmatrix") {
      axisval$z          = get(which_matrices[3])[,z_cols[1]]  
    }else{
      axisval$z    = apply(get(which_matrices[3])[,z_cols],1,mean)
    }
    return(axisval)
  })
  
  as.numeric.factor <- function(x) {match(x,levels(x))}
  if(is.factor(axisval$x)) axisval$x = as.numeric.factor(axisval$x)
  if(is.factor(axisval$y)) axisval$y = as.numeric.factor(axisval$y)
  if(is.factor(axisval$z)) axisval$z = as.numeric.factor(axisval$z)
  
  # Open 3d picture,
  open3d(...)
  
  # Get colours
  if(exists("obs.indv.colours",envir=forestFloor_graphics.env)) {
    colpal = get("obs.indv.colours",envir=forestFloor_graphics.env)
  } else {
    colpal = "black"
  }
  
  
  #should surfe of data also be plotted
  if(plot.surface) { 
    #compute grid around data
    get.seq = function(x) seq(min(x),max(x),length.out=grid.lines)
    XY = as.matrix(cbind(axisval$x,axisval$y),dimnames=NULL)
    ite.val=apply(XY,2,get.seq)
    gridXY=as.matrix(expand.grid(ite.val[,1],ite.val[,2]),dimnames=NULL) #grid coordinates
    g.points = grid.lines^2
    #     
    #     #rescale variables to equal to achieve equal influence in kNN-model
    #     sXY = scale(XY) #scale XY
    #     sgridXY = scale.by(scale.this=gridXY,by.this=sXY) #scale grid exactly as XY
    
    ##bootstrap knn estimated surface, giving gaussian-isch distance weigths
    out = kNN.surf(knnBag=knnBag,
                   XY=XY,
                   gridXY=gridXY,
                   k=k,
                   y=axisval$z,
                   bag.ratio = bag.ratio
    )
    
    #test the self-explainance of surface
    out2 = kNN.surf(knnBag=knnBag,
                    XY=XY,
                    gridXY=XY,
                    k=k,
                    y=axisval$z,
                    bag.ratio = bag.ratio
    )
    
    SSerrorSurf = sum((axisval$z-out2)^2)
    SSmodelSurf = sum((out2)^2)
    R2 = SSmodelSurf / (SSmodelSurf + SSerrorSurf)
    cat("fit of this surface plot is",round(R2,3),"% \n")
    
    # #debugging break
    # return(mget(ls()))
    
    xlab = names(x$X)
    
    
    #plot.surface
    persp3d(x=ite.val[,1], y=ite.val[,2], z=out,
            xlab = xyzlab$x,    ylab = xyzlab$y,    zlab = "feature contribution",
            aspect=c(1, 1, z_scale),main=paste("surf-fit=",round(R2,3),"%"),
            alpha=alpha.surf,col="#f2f2f2ff",
            ...)
    points3d(axisval$x,axisval$y,axisval$z,col=colpal,size=size.obs,alpha=alpha.obs,...)    
    
  }else{  # plot.surface = FASLSE, then data points only
    plot3d(axisval$x,axisval$y,axisval$z,
           col=colpal,
           aspect=c(1, 1, z_scale),
           size=size.obs,
           alpha=alpha.obs,
           xlab = xyzlab$x,
           ylab = xyzlab$y,
           zlab = "feature contribution",
           ...)
  }            
}




#f2 - show vec plot 2D and 3D
vec.plot = function(model,X,i.var,grid.lines=100,VEC.function=mean,zoom=1,limitY=F,forestFloor.col=FALSE) {

#compute grid range
d = length(i.var)
scales = lapply(i.var, function(i) {
  rXi = range(X[,i])
  span = abs(rXi[2]-rXi[1])*zoom/2
  center = mean(rXi)
  seq(center-span,center+span,length.out=grid.lines)
})

#expand grid range to a n-dimensional VEC-space and predict by model
anchor.points = as.matrix(expand.grid(scales),dimnames=NULL)
Xgeneralized=apply(X,2,VEC.function)    
Xtest.vec = data.frame(t(replicate(dim(anchor.points)[1],Xgeneralized)))
Xtest.vec[,i.var] = anchor.points
yhat.vec = predict(model,Xtest.vec)

#add observations to VEC-space
values.to.plot=X[,i.var]
Xmean=apply(X,2,VEC.function)
Xtest.obs = data.frame(t(replicate(dim(X)[1],Xgeneralized)))
Xtest.obs[,i.var] = values.to.plot
yhat.obs =  predict(model, Xtest.obs)

#defining one colour for all observations or import external colours
if(forestFloor.col) {col = with(forestFloor_graphics.env,obs.indv.colours)} else {col="#20202050"}

#plot VEC-space versus predictions (only 2D and 3D plot supported)
if(d==2) { #if 2D VEC space
  plot3d(x=values.to.plot[,1],y=values.to.plot[,2],z=yhat.obs,
         xlab=names(X)[i.var][1],ylab=names(X)[i.var][2],main="VEC-SURFACE",col=col)
  surface3d(x=scales[[1]],
            y=scales[[2]],
            z=yhat.vec,col="#404080",size=4,alpha=0.4)
}else{ #otherwise if 1D VEC-space
  if(limitY) {ylim = range(model$y)} else {ylim = NULL}
  plot(x=scales[[1]],y=yhat.vec,col="red",type="l",xlab=names(X)[i.var][1],ylim=ylim)
  points(values.to.plot,yhat.obs)
}
}



#f3 input a forestFloor object, computes convoluted feature contributions with kknn
#outout a new ff object as input with attached ff$FCfit
#kknn arguments can be accessed directly from this wrapper by userArgs.kknn
#if conflicting with wrappe
convolute_ff = function(ff,
                        these.vars=NULL,
                        k.fun=function() round(sqrt(n.obs)/2),
                        userArgs.kknn = alist(kernel="gaussian")
                        ) {
  
  n.obs=dim(ff$X)[1]
  n.vars=dim(ff$X)[2]
  k=k.fun()
  if(is.null(these.vars)) these.vars = 1:n.vars
  
  #merge user and wrapper args
  defaultArgs.kknn = alist(formula=fc~.,data=Data,kmax=k,kernel="gaussian")
  kknn.args=append.overwrite.alists(userArgs.kknn,defaultArgs.kknn)
  
  #iterate for seleceted variables
  ff$FCfit = sapply(these.vars, function(this.var) {
    Data = data.frame(fc=ff$FCmatrix[,this.var],x=ff$X[,this.var])
    knn.obj = do.call("train.kknn",kknn.args)$fitted.values
    knn.obj[[length(knn.obj)]]
  })
  ff  
}


#f3 input a forestFloor object, as convolute_ff but do not iterate all variables. Instead wrapper will use
#kknn to convolute a designated featureContribution with designated features - oftenly the corresponding features.
convolute_ff2 = function(ff,
                        Xvars,
                        FCvars = NULL,
                        k.fun=function() round(sqrt(n.obs)/2),
                        userArgs.kknn = alist(kernel="gaussian")
) {
  if(is.null(FCvars)) FCvars = Xvars
  n.obs=dim(ff$X)[1]
  k=k.fun()
    
  #merge user and wrapper args
  defaultArgs.kknn = alist(formula=fc~.,data=Data,kmax=k,kernel="gaussian")
  kknn.args=append.overwrite.alists(userArgs.kknn,defaultArgs.kknn)
  
  #collect coloumns
  if(length(FCvars)>1) {
    fc = apply(ff$FCmatrix[,FCvars],1,sum)
  } else {
    fc = ff$FCmatrix[,FCvars]
  }
  x = ff$X[,Xvars]
  
  #compute topology
    Data = data.frame(fc=fc,x=x)
    knn.obj = do.call("train.kknn",kknn.args)$fitted.values
    out = knn.obj[[length(knn.obj)]]
}

#f4 input a forestFloor object, as convolute_ff but do not iterate all variables. Instead wrapper will use
#kknn to convolute a designated featureContribution with designated features - oftenly the corresponding features.
convolute_grid = function(ff,
                          Xvars,
                          FCvars = NULL,
                          grid = 30,
                          limit = 3,
                          zoom = 3,
                          k.fun=function() round(sqrt(n.obs)/2),
                          userArgs.kknn = alist(kernel="gaussian")
) {
  
  #input defaults
  if(is.null(FCvars)) FCvars = Xvars
  n.obs=dim(ff$X)[1]
  k=k.fun() # will be overided if a "k=" argument is provided in userArgs.kknn-list
  
  #collect coloumns of data
  if(length(FCvars)>1) {
    fc = apply(ff$FCmatrix[,FCvars],1,sum)
  } else {
    fc = ff$FCmatrix[,FCvars]
  }
  X = ff$X[,Xvars]
  
  #make grid, or use external grid if provided
  if(length(grid)==1) {
    X = box.outliers(X,limit=limit,normalize=F)
    get.seq = function(x) {
      upper = (max(x) - mean(x)) * zoom
      lower = (min(x) - mean(x)) * zoom
      seq(lower+mean(x), upper+mean(x),length.out=grid)
    }
    ite.val=lapply(1:dim(X)[2],function(i) get.seq(X[,i])) #for each variable define span of grid
    gridX = as.data.frame(expand.grid(ite.val))
    names(gridX) = names(X)
    
    ##WOUPS NEVERMIND it worked...
#     #could not make list of vectors work with expand.grid as doc promises
#     #this hack calls expand.grid with specified arguments
#     
#     #hack - make charecter string of appropiate code
#     run.this = "alist("
#     for(i in 1:dim(X)[2]) run.this = paste(run.this,names(X)[i]," = ite.val[,",i,"],",sep="")
#     run.this = substr(run.this,1,nchar(run.this)-1)
#     run.this = paste(run.this,")")
#     #hack - parse and evaluate string into a argument list of one argument for each dimension
#     arg.list = eval(parse(text=run.this))
#     #call expand.grid with specified arguments
#     gridX=as.data.frame(do.call(expand.grid,arg.list)) #grid coordinates

  } else {
    gridX = grid
  }
  
  #prepare data
  Data = data.frame(fc=fc,x=ff$X[,Xvars])
  gridX = data.frame(x=gridX)
  
  #merge user args and args of this wrapper function, user args have priority
  defaultArgs.kknn = alist(formula=fc~.,train=Data,k=k,kernel="gaussian",test=gridX)
  kknn.args=append.overwrite.alists(userArgs.kknn,defaultArgs.kknn)
    
  #execute kknn function and retrive convolution
  gridFC = do.call("kknn",kknn.args)$fitted.values
  
  #gather results in data frame, with convoluted feature contributinos in grid
  return(data.frame(fc = gridFC,gridX))
}



#sf5 scale data and grid, to allow knn 
scale.by = function(scale.this,by.this) {
  center = attributes(by.this)$'scaled:center'
  scales = attributes(by.this)$'scaled:scale'
  nvars = dim(scale.this)[2]
  sapply(1:nvars, function(i) (scale.this[,i]-center[i])/scales[i])
}




#sf6  reduce outliers to within limit of 1.5 std.dev and/or output as normalized 
box.outliers = function(x,limit=1.5,normalize=T) {
  
  sx=scale(x)
  if(limit!=FALSE) {
    sx[ sx>limit] =  limit
    sx[-sx>limit] = -limit
  }
  
  if(normalize) { 
    sx.span = max(sx) - min(sx)
    sx = sx - min(sx)
    sx = sx / sx.span
  } else {
    obs=attributes(sx)$"dim"[1]
    if(dim(sx)[2]>1) {
    sx = sx * t(replicate(obs,attributes(sx)$"scaled:scale")) +
              t(replicate(obs,attributes(sx)$"scaled:center"))
    } else {
    sx  = sx * attributes(sx)$"scaled:scale" + attributes(sx)$"scaled:center" 
    }
  }

  if(class(x)=="data.frame") {
    sx = as.data.frame(sx,row.names=row.names(x))
    names(sx) = names(x)
  }
  return(sx)
}

# #sf7  - calculate deviance from surface
# surf.error = function() {
#   outs=replicate(knnBag, {  #the following is replicated/performed severeal ~20 times
#     this.boot.ind = sample(dim(XY)[1]*bag.ratio,replace=T) #pick a bootstrap from samples             
#     sXY = scale(XY[this.boot.ind,])  #scale bootstrap to uni-variance
#     sgridXY = scale.by(scale.this=gridXY,by.this=sXY) #let grid be scaled as this bootstrap was scaled to sXY
#     out=knn.reg(train=sXY,
#                 test=sgridXY,
#                 y=axisval$z[this.boot.ind],
#                 k=k,
#                 algorithm="kd_tree")$pred  #predict grid from bootstrap of samples
#   })
# }

#sf7 estimate surface with kNNbag, depends on "sf5 - scale.by"
kNN.surf = function(knnBag,
                    XY,
                    gridXY,
                    k,
                    y,
                    bag.ratio=.8,
                    replace=T) {
  outs=replicate(knnBag, {  #the following is replicated/performed severeal ~20 times
    this.boot.ind = sample(dim(XY)[1]*bag.ratio,replace=T) #pick a bootstrap from samples             
    sXY = scale(XY[this.boot.ind,])  #scale bootstrap to uni-variance
    sgridXY = scale.by(scale.this=gridXY,by.this=sXY) #let grid be scaled as this bootstrap was scaled to sXY
    out=knn.reg(train=sXY,
                test=sgridXY,
                y=y[this.boot.ind],
                k=k,
                algorithm="kd_tree")$pred  #predict grid from bootstrap of samples
  })
  out = apply(outs,1,mean) # collect predictions
  return(out)
}



#sf8 neat function to help increase adaptability of wrappers, default args defined by wrapper. User can 
#input an alist and by this function the wrapper will append new args and overwrite conflicting arguments.
#one set of args either user or defaults are set as master
append.overwrite.alists= function(masterArgs,slaveArgs) {
  slaveArgs.to.overwrite = names(slaveArgs) %in% names(masterArgs)
  for(i in which(slaveArgs.to.overwrite))  slaveArgs[i] = masterArgs[match(names(slaveArgs[i]),names(masterArgs))]
  masterArgs.to.append = !(names(masterArgs) %in% names(slaveArgs))
  c(slaveArgs,masterArgs[masterArgs.to.append])
}



forestFloor = function(rfo,X,calc_np=FALSE) { 
  
  #check the RFobject have a inbag
  if(is.null(rfo$inbag)) stop("input randomForest-object have no inbag, set keep.inbag=T,
try, randomForest(X,Y,keep.inbag=T) for regression where Y is numeric
and, cinbag(X,Y,keep.inbag=T,keep.forest=T) for binary-class where Y is factor
..cinbag is from trimTrees package...
error condition: if(is.null(rfo$inbag))")
   
  #make node status a integer matrix
  ns = rfo$forest$nodestatus
  storage.mode(ns) = "integer"
   
  
  #translate binary classification RF-object, to regression mode
  if(rfo$type=="classification") {
    if(length(levels(rfo$y))!=2) stop("no multiclass, must be binary classification.
                                      error condition: if(length(levels(rfo$y))!=2")
    print("RF is classification, converting factors/categories to numeric 0 an 1")
    Y = as.numeric((rfo$y))-1
    cat(" defining",levels(rfo$y)[1]," as 0\n defining",levels(rfo$y)[2],"as 1")
    rfo$forest$leftDaughter  = rfo$forest$treemap[,1,] #translate daughter representation to regression mode
    rfo$forest$rightDaughter = rfo$forest$treemap[,2,] 
    ns[ns==1] = -3  ##translate nodestatus representation to regression mode
    if(is.null("rfo$inbagCount")) stop("classification topology not supported with randomForest() {randomForest}
Grow forest with cinbag::trimTrees instead of randomForest(). The two
functions are identical, except cinbag() entails a more detailed inbag record,
which is needed to estimate binary node probabilities.
error condition:  if(is.null('rfo$inbagCount'))")
    
    if(!calc_np) stop("node predictions must be re-calculated for random forest of type classification, set calc_np=T)
error conditions: if(!calc_np && rfo$type='classification')")
    
    inbag = rfo$inbagCount
    } else {
    Y=rfo$y
    inbag = rfo$inbag
  }


  #preparing data, indice-correction could be moved to C++
  #a - This should be fethed from RF-object, flat interface
  ld = rfo$forest$leftDaughter-1 #indice correction, first element is 0 in C++ and 1 in R.
  storage.mode(ld) = "integer"
  rd = rfo$forest$rightDaughter-1
  storage.mode(rd) = "integer"
  bv = rfo$forest$bestvar-1
  storage.mode(bv) = "integer"
  np = rfo$forest$nodepred
  storage.mode(np) = "double"
  bs = rfo$forest$xbestsplit
  storage.mode(bs) = "double"
  ib = inbag
  storage.mode(ib) = "integer"
  Yd = as.numeric(Y)
  storage.mode(Yd) = "double"
  ot  = rfo$oob.times
  storage.mode(ot) = "integer"
 
 
  ##recording types of variables
  xlevels = unlist(lapply(rfo$forest$xlevels,length),use.names=F)
  xl = xlevels
  storage.mode(xl) = "integer"
  varsToBeConverted = xlevels>1

  ##Converting X to Xd, all factors change to level numbers
  Xd=X
  for(i in 1:dim(Xd)[2]) {
    if(varsToBeConverted[i]) {
      Xd[,i] = as.numeric(Xd[,i])-1  
    }
  }  
  Xd=as.matrix(Xd)
  storage.mode(Xd) = "double"
  
  #outout variable
  localIncrements = Xd*0
  storage.mode(localIncrements) = "double"
  
  #should activities of nodes be reestimated(true) or reused from randomForest object(false)
  calculate_node_pred=calc_np
    
  # C++ function, recursively finding increments of all nodes of all trees
  # where OOB samples are present. vars, obs and ntree is "passed by number"
  # Anything else is passed by reference. Found increments are imediately
  # summed to localIncrements matrix.
  recTree(
    #passed by number
    vars=dim(X)[2], 
    obs=dim(X)[1],             
    ntree=rfo$ntree,
    calculate_node_pred=calculate_node_pred,
    #passed by reference
    X=Xd,  #training data, double matrix [obs,vars] 
    Y=Yd,
    leftDaughter = ld,  #row indices of left subnodes, integer matrix [nrnodes,ntree] 
    rightDaughter = rd, #...
    nodestatus = ns,    #weather node is terminal or not,      
    xbestsplit = bs,          
    nodepred = np,          
    bestvar = bv,
    inbag = ib,
    varLevels = xl,
    ot,  #oob.times
    localIncrements = localIncrements #output is written directly to localIncrements from C++
  )
  
  
  
#writing out list
  imp = as.matrix(rfo$importance)[,1]
  out = list(X=X,Y=Y,
             importance = imp,
             imp_ind = sort(imp,decreasing=T,index.return=T)$ix,
             FCmatrix = localIncrements
  )
  class(out) = "forestFloor"
  return(out)
}