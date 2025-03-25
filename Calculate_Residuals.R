Calculate_Residuals<-function(y,x,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred){ #A function that calculates the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation.
  if (j==1){
    indexes=Cell_Indexes(x,Tess[[j]],Dim[[j]]);
    CurrentTessPred<-Pred[[j]][indexes]
    SumOfAllTess=SumOfAllTess-CurrentTessPred}
  else{
    indexes=Cell_Indexes(x,Tess[[j]],Dim[[j]]);
    CurrentTessPred<-Pred[[j]][indexes]
    SumOfAllTess=SumOfAllTess+LastTessPred-CurrentTessPred;
  }
  
  IndexesStar=Cell_Indexes(x,TessStar[[j]],DimStar[[j]]);
  R_j<-y-SumOfAllTess
  
  #Initializing Sizes
  
  R_ijOld=rep(0,length(Pred[[j]]))
  n_ijOld=rep(0,length(Pred[[j]]))
  
  for (i in 1:length(Pred[[j]])){
    R_ijOld[i]<-sum(R_j[indexes==i])
    n_ijOld[i]<-sum(indexes==i)
  }
  
  R_ijNew=rep(0,length(TessStar[[j]][,1]))
  n_ijNew=rep(0,length(TessStar[[j]][,1]))
  
  for (i in 1:length(TessStar[[j]][,1])){
    R_ijNew[i]<-sum(R_j[IndexesStar==i])
    n_ijNew[i]<-sum(IndexesStar==i)}
  
  return(list(R_ijOld,n_ijOld,R_ijNew,n_ijNew,SumOfAllTess,IndexesStar,indexes))
}
