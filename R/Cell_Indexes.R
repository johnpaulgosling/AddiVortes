library('FNN')

Cell_Indexes<-function(x,Tess,Dim){ #Gives the row (the center) of the tessellation that each obseravtion falls within.
  if (length(Tess[,1])==1){ #only 1 centre
    CellsForGivenTess=rep(1,length(x[,1]))
  }
  else{ #multiple
    CellsForGivenTess=knnx.index(Tess,matrix(x[,Dim],ncol = length(Dim)),1)
  }
  return(CellsForGivenTess)
}
