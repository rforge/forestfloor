#load libaries

load_stuff = function() {
  library(rfFC)
  library(foreach)
  library(randomForest)
  library(doSNOW)
}


### function to force close a lost cluster, where keys were saved globally
globalCloseCluster = function() stopCluster(get("reg.cl"))

###function to compute randomForest cross validated predictions and feature contributions  

##methods

#m1 print output
print.forestFloor = function(RFCV_object) {
  print("RFCV_contribution model:")
  print(c("CV result of repititions"))
  print(paste("mean validation Q² =",round(mean(RFCV_object$Q2_reps),2)))
  print(paste("sd of validation Q² =",round(sd(RFCV_object$Q2_reps),6)))
  
  print(paste("mse validated =",round(mean(RFCV_object$mse_reps),2)))
  print(paste("mad validated =",round(mean(RFCV_object$mad_reps),2)))
  dims=dim(RFCV_object$predcube)
  print(paste("obs =",dims[1]," ,fold =",dims[2]," ,reps = ",dims[3]))
  print("mean importance from cross-validations, increase MSE%")
  print(round(data.frame(RFCV_object$importance[RFCV_object$imp_ind]),2))
}

#m2 plot output
plot.forestFloor = function(RFCV_object,
                            plot_seq=NULL,
                            colour_by=1,
                            col_axis = 2,
                            alpha="auto"
                            )
  {
  
  
  
  
  #save previous pars
  pars = par(no.readonly = TRUE)
  par(mar=c(2,2,1,1))
  
  
  #short for phys.val and feature contribution in object
  X = RFCV_object$X
  FCs = RFCV_object$FCsCV
  
  
  #setting transparancy variable
  if(alpha=="auto") alpha = min(max(400/dim(X)[1],0.2),1)
  
  #choosing to plot first 18
  if(is.null(plot_seq)) plot_seq = 1:min(dim(X)[2],24)
  
  #make everything numeric, save jitter.template
  jitter.template=rep(FALSE,dim(X)[2])
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
  imp = RFCV_object$importance     #fetch importance
  imp.ind = RFCV_object$imp_ind    #fetch importance ranking/indices
  
  
  
  #set default.colour
  colours = "black"
  
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
  
  if(colour_by=="black") colours="black"
  
  #Save this colouring globally, for later 3D plotting
  assign("global.col",colours,.GlobalEnv)
  
  
  ##plot the n.plots most important variables
  
  for(i in plot_seq) {
    plot(
      data.frame(
        physical.value        = jitter(X[,imp.ind[i]],factor=jitter.template[imp.ind[i]]*2),
        partial.contribution  = FCs[,imp.ind[i]]),
      main = names(imp)[imp.ind[i]],
      cex= 0.5,
      col = colours)
  }
  
  par(pars)
}

#m3 3d show function
show3d = function(RFCV_object,
         order_by_importance=F,
         which_matrices=c("X","X","FCsCV"),
         x_cols=2,y_cols=11,z_cols=c(2),
         plot.surface=T,
         grid.lines=30,
         k=5,
         alpha=.4,
         z_scale=.7,
         knnBag=20,
         bag.ratio=0.5) {
  library(rgl)
  
  #apply funciton later requires two inputs
 
  
  cols.fix = function(i)if(length(i)<2){i=c(i,i)}else{i}
  x_cols = cols.fix(x_cols)
  y_cols = cols.fix(y_cols)
  z_cols = cols.fix(z_cols)
  
  
  axisval = with(RFCV_object,{
    axisval=c()
    if(order_by_importance) {
      imp.ind = sort(importance,decreasing=T,index.return=T)$ix
      ind_by_imp = function(cols,imp.ind) imp.ind[cols]
      x_cols = ind_by_imp(x_cols,imp.ind)
      y_cols = ind_by_imp(y_cols,imp.ind)
      z_cols = ind_by_imp(z_cols,imp.ind)
    }
      
      if(which_matrices[1]=="X") {
        axisval$x = get(which_matrices[1])[,x_cols[1]]
      }else{
        axisval$x = apply(get(which_matrices[1])[,x_cols],1,mean)
      }
      if(which_matrices[2]=="X") {
        axisval$y = get(which_matrices[2])[,y_cols[1]]
      }else{
        axisval$y = apply(get(which_matrices[2])[,y_cols],1,mean)
      }
      if(which_matrices[3]=="X") {
        axisval$z = get(which_matrices[3])[,z_cols[1]]
      }else{
        axisval$z = apply(get(which_matrices[3])[,z_cols],1,mean)
      }
    
    return(axisval)
  })
  
  as.numeric.factor <- function(x) {match(x,levels(x))}
  if(is.factor(axisval$x)) axisval$x = as.numeric.factor(axisval$x)
  if(is.factor(axisval$y)) axisval$y = as.numeric.factor(axisval$y)
  if(is.factor(axisval$z)) axisval$z = as.numeric.factor(axisval$z)

  print(axisval$x)
  print(axisval$y)
  print(axisval$z)
  
  # Open 3d picture,
  open3d()
  bg3d("white")
  material3d(col="black")  
  if(exists("global.col",env=.GlobalEnv)) {
    colpal = get("global.col")
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
    
    #rescale variables to equal to achieve equal influence in kNN-model
    sXY = scale(XY) #scale XY
    sgridXY = scale.by(scale.this=gridXY,by.this=sXY) #scale grid exactly as XY
      
    ##bootstrap knn estimated surface, giving gaussian-isch distance weigths
    require(FNN)
    outs=replicate(knnBag, {  
      sXY = scale(XY)
      sgridXY = scale.by(scale.this=gridXY,by.this=sXY)
      this.boot.ind = sample(dim(XY)[1]*bag.ratio,replace=T)             
      out=knn.reg(train=sXY[this.boot.ind,],
                    test=sgridXY,
                    y=axisval$z[this.boot.ind],
                    k=k,
                    algorithm="kd_tree")$pred
      })
      out = apply(outs,1,mean) # collect predictions
      
      #plot.surface
      persp3d(x=ite.val[,1], y=ite.val[,2], z=out,
              xlab = "X",    ylab = "Y",    zlab = "partial functions",
              aspect=c(1, 1, z_scale),
              alpha=alpha,col="#f2f2f2ff")
      points3d(axisval$x,axisval$y,axisval$z,col=colpal)    
  
  }else{  # plot.surface = FASLSE, then data points only
    points3d(axisval$x,axisval$y,axisval$z,col=colpal,aspect=c(1, 1, z_scale))
  }            
}



#sf2
mad = function(x,y) mean(abs(x-y))
#sf3
mse = function(x,y) mean((x-y)^2)
#sf4 classification
acc = function(x,y) mean(x==y)

#sf5 scale data and grid, to allow knn 
scale.by = function(scale.this,by.this) {
  center = attributes(by.this)$'scaled:center'
  scales = attributes(by.this)$'scaled:scale'
  nvars = dim(scale.this)[2]
  sapply(1:nvars, function(i) (scale.this[,i]-center[i])/scales[i])
}




#sf6  reduce outliers to within limit of 1.5 and normalize
box.outliers = function(x,limit=1.5) {
  x=scale(x)
  x[ x>limit] =  limit
  x[-x>limit] = -limit
  #normalize component to [0;1]
  x=x-min(x)
  x=x/ (limit*2)
  return(x)
}



##function
forestFloor = function(X,Y,
                             .ntree  = 150,
                             fold   =   10,
                             reps   =   5,
                             parallel = T,
                             n_jobs =  5,
                             spawn_jobs_by_reps = F,
                             single_model_parallel = F,
                             classwt = NULL,
                             obs_maxnodes_ratio = 1,
                             replace=F
) {
  
  ###function to compute randomForest cross validated predictions and feature contributions  
  #function supports parallel processing
  
  print(reps)
  
  ## Setting intructions for parallel computing and handling export of packages, if needed
  if(parallel) {
    if (n_jobs==-1) n_jobs = detectCores()
    cl = makeCluster(n_jobs)
    assign("reg.cl",cl,env=.GlobalEnv)# in case function chrashed, ' the keys' for cluster is saved in global environment
    registerDoSNOW(cl)
    if(spawn_jobs_by_reps) {
      #parallel by reps, export packages in repitition loop
      '%doReps%' = get("%dopar%")
      '%doFold%' = get("%do%")
      reps_packages = c("randomForest","rfFC","foreach","trimTrees") #packages exported to cores
      fold_packages = NULL #if not parallel, the package is still in scope, no export needed
      reps_export=NULL
    } else {
      #parallel by fold, export used packages in fold loop
      '%doReps%' = get("%do%")
      '%doFold%' = get("%dopar%")
      reps_packages = NULL #if not parallel, the package is still in scope, no export needed
      fold_packages = c("randomForest","rfFC","trimTrees") #packages to exported to cores, foreach not needed
    }
  } else {
    # if not parallel, do everything sequential, all packages are already in scope
    '%doReps%' = get("%do%")
    '%doFold%' = get("%do%")
    reps_packages = NULL
    fold_packages = NULL
  }
  
  ##Here start reps repitition of validation loops, .packages pass designated libraries to cores
  CV_result = foreach(this.rep = 1:reps,.combine=c,.packages = reps_packages) %doReps% {
                        
    #Print progress
    print(paste(c("commencing repitition:",this.rep)))
    
    ##make new partion keys for this n-Fold cross validation repitition
    this.rep.shuf     = sample(1:dim(X)[1],dim(X)[1])      # make a shuffle key, 0 3 1 6 2 8 4 7 5 9 ...
    this.rep.part     = (0:(dim(X)[1]-1)%%fold)+1          # make a partion key, 1 2 3 4 1 2 3 4 1 2 ...
                        
    ##Here starts nfold validation loop, .combine ensures all result objects are combined into vector, "c"
    this.rfFCCV = foreach(this.fold = 1:fold,.combine=c,.packages=fold_packages) %doFold% {
                                                
      ##define training- and testsets
      test = this.rep.shuf[this.rep.part==this.fold]  #mix keys for each fold. This yields test indices
      this.Ytrain = Y[-test]
      this.Xtrain = X[-test,]
      this.Xtest  = X [test,]
      this.Ytest  = Y [test]
                                            
      ##train 1 randomForest model for one fold for one repitition
      .maxnodes  = ceiling(length(this.Ytrain)/obs_maxnodes_ratio)
      
      
      #using modifoed randomForest to caluculate inbagCount if replace == TRUE
      if(replace==T) {
        func.rf = cinbag
      } else {
        func.rf = randomForest
      }
      
      RF_model = func.rf(this.Xtrain,
                              this.Ytrain,
                              keep.inbag = T,
                              classwt    = classwt,
                              replace    = replace,
                              importance = T, 
                              ntree      = .ntree,
                              maxnodes   = .maxnodes)
                            
      ##compute forest contributions and insert in zero_matrix with dimensions as X
      #these two major functionsis suppurted by the rfFC-package
      RF_model$inbag = RF_model$inbagCount
      RF_increments = getLocalIncrements(RF_model,this.Xtrain)
      RF_contributions = featureContributions(RF_model,RF_increments,this.Xtest)
      zero_matrix = matrix(0, length(X),nrow=dim(X)[1])
      zero_matrix[test,]  = RF_contributions # put results in to matrix
      RF_contributions_matrix = zero_matrix
                            
      ##this fold this rep prediction for cross validation
      this.Ypred = predict(RF_model,this.Xtest) #predict testSet for cross validation purpose
      #combine prediction and true value and observation-,fold- and reps- number into prediction table
      predtable = data.frame(
        obs.ind = test,
        pred = this.Ypred,
        true = this.Ytest,
        this.rep  = rep(this.rep,length(this.Ypred)),
        this.fold = rep(this.fold,length(this.Ypred))
      )
                            
      #return combined contributions, predtable and measured importance as objects in vector
      out = c(list(RF_contributions_matrix),
               list(predtable),
               list(RF_model$importance[,1]))
      
#       class_list = c("RF_contributions_matrix",
#                      "predtable",
#                      "RF_model$importance")
#       
#       #assigning classes to outputs, for later reconstruction
#       for(i in 1:length(out)) class(out[[i]])=class_list[i]
      return(out)
      
               #list(RF_model_list
      } # here stops fold loop
  } #here stops reps loop

  #build one non cross-validated model
  if(parallel) {
    '%run%'=get("%dopar%")
  }else{
    '%run%'=get("%do%")
  }
  print("finished cross validation...")
  
  .maxnodes = ceiling(length(Y) / obs_maxnodes_ratio)
  
  if(single_model_parallel) {
    full_RF_model = foreach(sub_tree = rep(ceiling(.ntree/5),n_jobs),
                            .combine = combine,
                            .packages ="randomForest",
                            importance=T) %run% {
      randomForest(X,Y,ntree=sub_tree,maxnodes=.maxnodes)
    }
  } else {
    full_RF_model = randomForest(X,Y,
                                 ntree=.ntree,
                                 maxnodes=.maxnodes,
                                 importance=T)
  }
  
  
  if(parallel) stopCluster(cl)

 print("closing cluster")
  
  ##Collecting and rearranging results from jobs
  #calculate indices to separate objects in very long result vector
  all_indices = 1:(fold*reps*3)
  rfFCCVs_ind = which(all_indices%%3==1)
  predtable_ind = which(all_indices%%3==2)
  importance_ind = which(all_indices%%3==0)
  
#   #get models
#   print("so far!!")
#   RF_model_list = CV_result[RF_model_ind] # get interlaced results, of RF_models
#   combined_RF_model = foreach(i = 1:(fold*reps), .combine=combine) %do% {
#     unlisted_model =  RF_model_list[[1]]
#   }
  


  #simplify FCs to array of [obs,vars,reps]
  print("combine feature contrinutions")
  rfFCCVs_list = CV_result[rfFCCVs_ind] # get interlaced results, of rfFC
  rfFCCVs_array = simplify2array(rfFCCVs_list) #convert to array of dim [obs,vars,(repsXfold)]
  FCsCV = apply(rfFCCVs_array,1:2,sum)/reps #combine folds and get mean FCs across repetitions 

  #simplify importance to mean importance
  print("agregating importance")
  importance_list = CV_result[importance_ind]# get interlaced results, of importance 
  importance_array = simplify2array(importance_list)#
  importance_mean = apply(importance_array,1,mean)
  importance_indices = sort(importance_mean,decreasing=T,index.return=T)$ix
  #plot(importance_mean[importance_indices])
  
  #concatenate predtable to large table
  print("combining predtables")
  predtable_list = CV_result[predtable_ind]# get interlaced results, of importance 
  predtable_df = do.call(rbind,predtable_list)
  predtable_replist = lapply(1:reps,function(this.rep) {
    predtable_rep=predtable_df[this.rep==predtable_df$this.rep,]
    out = predtable_rep[sort(predtable_rep[,1],index.return=T)$ix,]
    simplify2array(out)
  })
  predtable_array = simplify2array(predtable_replist)
  
  
  #produce CV statistics
  print("computing cv statistics")
  do.table2 =  function(func,col1,col2) { #function to iterate pred.table
    lapply(1:reps,{
      function(this.rep) {
        do.call(
          what=func,
          args = list(predtable_array[,col1,this.rep],predtable_array[,col2,this.rep])
        )
      }
    })
  }
  
  Q2_reps  = unlist(do.table2(cor,2,3))^2
  mse_reps = unlist(do.table2(mse,2,3))
  mad_reps = unlist(do.table2(mad,2,3))
  

  #writing out structure
  print("wrapping up...")
  out = list(X=X,Y=Y,
             predcube = predtable_array,
             importance = importance_mean,
             imp_ind = importance_indices,
             FCsCV = FCsCV,
             mse_reps = mse_reps,
             Q2_reps  = Q2_reps,
             mad_reps = mad_reps,
             one_RF_model = full_RF_model
             )

  class(out) = "forestFloor"
  return(out)
}



#get surf - standalone plotting of surface and points
kNN.surf = function(x,y,z,
                    grid.lines=30,
                    k=5,
                    func=mean,
                    alpha=.4,
                    z_scale=.1,
                    knnBag=20,
                    bag.ratio=0.5,
                    col=NULL){#= col="global" will fetch global colours  
  get.seq = function(x) seq(min(x),max(x),length.out=grid.lines)
  XY = as.matrix(cbind(x,y),dimnames=NULL)
  ite.val=apply(XY,2,get.seq)
  gridXY=as.matrix(expand.grid(ite.val[,1],ite.val[,2]),dimnames=NULL)
  g.points = grid.lines^2
  
  # scale variables
  scale.by = function(scale.this,by.this) {
    center = attributes(by.this)$'scaled:center'
    scales = attributes(by.this)$'scaled:scale'
    nvars = dim(scale.this)[2]
    sapply(1:nvars, function(i) (scale.this[,i]-center[i])/scales[i])
  }
  sXY = scale(XY) #scale XY
  sgridXY = scale.by(scale.this=gridXY,by.this=sXY) #scale grid as XY
  
  ##bootstrap knn estimated surface, giving gaussian-isch distance weigths
  require(FNN)
  outs=replicate(knnBag, {  
    sXY = scale(XY)
    sgridXY = scale.by(scale.this=gridXY,by.this=sXY)
    this.boot.ind = sample(dim(XY)[1]*bag.ratio,replace=T)             
    out=knn.reg(train=sXY[this.boot.ind,],
                test=sgridXY,
                y=z[this.boot.ind],
                k=k,
                algorithm="kd_tree")$pred
  })
  out = apply(outs,1,mean)
  print(col)
  if((col=="global" && exists("global.col",env=.GlobalEnv))) {
      colpal = get("global.col")
    } else {
      colpal = col
      if(is.null(col)) colpal = "black"
    }
    
    
  open3d()
  bg3d("white")
  material3d(col="black")
  persp3d(ite.val[,1], ite.val[,2], out, aspect=c(1, 1, z_scale),
          xlab = "X", ylab = "Y", zlab = "partial function",alpha=alpha,col="#ffc2c2")
  points3d(XY[,1],XY[,2],z,col=colpal)
}  

