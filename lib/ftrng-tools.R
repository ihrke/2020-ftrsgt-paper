evans.rng <- function(x, nitems){
  T=length(x)
  M=matrix(0, nrow=nitems,ncol=nitems)
  for(i in 1:nitems){
    for(j in 1:nitems){
      for(t in 1:(T-1)){
        if(x[t]==i & x[t+1]==j){
          M[i,j]=M[i,j]+1
        }
      }
    }
  }
  y=as.vector(M)
  num=sum(y*log(y), na.rm=T)
  y=apply(M, 1, sum)
  denom=sum(y*log(y), na.rm=T)
  num/denom
}
