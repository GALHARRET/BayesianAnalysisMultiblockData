algoMCOAv=function(X,R=2){
  Xt=list()
  K<-length(X)
  for(k in 1:K) Xt[[k]]=as.data.frame(X[[k]])
  XX=ade4::ktab.list.df(Xt)
  res=ade4::mcoa(XX, scannf = FALSE, nf = R,option = "uniform")
  U=as.matrix(res$SynVar)
  V<-list()
  SIG<-list()
  for(k in 1:K){
    V[[k]]=t(X[[k]])%*%U
    SIG[[k]]=sqrt(diag(t(U)%*%X[[k]]%*%V[[k]]))
  }
  I=sort(apply(matrix(unlist(SIG)^2,ncol=K),1,mean),decreasing = T,index.return=T)$ix
  U=U[,I]
  for(k in 1:K){
    SIG[[k]]=SIG[[k]][I]
    V[[k]]=V[[k]][,I]
  }
  var_explained=matrix(NA,nrow=K,ncol=R)
  rownames(var_explained)=paste("Block",1:K,sep="")
  colnames(var_explained)=paste("Dim",1:R,sep="")
  for(k in 1:K) var_explained[k,]=100*SIG[[k]]^2/sum(SIG[[k]]^2)
  return(list(U=as.matrix(res$SynVar),SIG=SIG,V=V,var_explained=var_explained))
}


