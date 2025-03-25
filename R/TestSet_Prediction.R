TestSet_Prediction<-function(x,m,Tess,Dim,Pred){ #A function that derives a prediction for the output variable give an observation, for one iteration in AddiVortes.
  Prediction=rep(0,length(x[,1]));
  for (j in 1:m){
    NewTessIndexes=Cell_Indexes(x,Tess[[j]],Dim[[j]]);
    Prediction=Prediction+Pred[[j]][NewTessIndexes]
  }
  return(Prediction)
}
