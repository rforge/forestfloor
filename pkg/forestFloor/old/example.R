#simulate data
obs=3000 
vars = 6 

X = data.frame(replicate(vars,rnorm(obs)))
Y = with(X, X1^2 + sin(X2*pi) + 2 * X3 * X4 + 0.2 * rnorm(obs))


#grow a forest, remeber to include inbag
rfo=randomForest(X,Y,keep.inbag = TRUE)


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
show3d(ff,x_cols=3,y_cols=4,z_cols=3:4,k=15) 


#although no interaction, a joined additive effect of X1 and X2
#can also be informative to display in 3D
plot(ff,plot_seq=NULL,colour_by=4,col_axis=2) #use plot first to define colours \cr
show3d(ff,x_cols=1,y_cols=2,z_cols=1:2) 
