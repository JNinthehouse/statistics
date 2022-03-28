mu.test.multinorm=function(x,y=NULL,mu0=rep(0,ncol(x)),Sigma0=-1){
  # single population : Sigma0 is known and Sigma0 is unknown
  # two populations   : Sigma0 is unknown 
  
  if (is.null(y)){
    x=as.matrix(x)
    n=nrow(x)
    p=ncol(x)
    
    if (Sigma0!=-1){
      X.bar=apply(x, 2, mean)
      T1=n*t(X.bar-mu0)%*%solve(Sigma0)%*%(X.bar-mu0)
      
      a2=qchisq(1-alpha, p)
      
      reject=matrix(c(T1, a2), nrow=1)
      rownames(reject)=c("Reject")
      colnames(reject)=c("Obs", ">1-alpha")
      
      pv=round(1-pchisq(T1, p),3)
      return(list(Reject.area=reject, Xmean=X.bar, Chi.obs=T1, p.value=pv))}
    
    else if (Sigma0==-1){
      X.bar=apply(x, 2, mean)
      A=(n-1)*var(x)
      
      T2=(n-1)*n*t(X.bar-mu0)%*%solve(A)%*%(X.bar-mu0)
      FF=(n-p)/((n-1)*p)*T2
      
      pv=round(1-pf(FF, p, n-p),3)
      return(list(Xmean=X.bar, F.obs=FF, p.value=pv))
    }
  }
  
  
  else if (!is.null(y)){
    x=as.matrix(x)
    y=as.matrix(y)
    n1=nrow(x)
    n2=nrow(y)
    p=ncol(x)
    if (n1==n2){
      X.bar=apply(x, 2, mean) 
      A1=(n1-1)*var(x)
      Y.bar=apply(y, 2, mean)
      A2=(n2-1)*var(y) 
      A=(A1+A2)/(n1+n2-2)
      
      T2=(n1*n2/(n1+n2))*t(X.bar-Y.bar)%*%solve(A)%*%(X.bar-Y.bar)
      FF=(n1+n2-2-p+1)/((n1+n2-2)*p)*T2
      
      pv=round(1-pf(FF, p, (n1+n2-p-1)),3)
      return(list(Xmean=X.bar, Ymean=Y.bar, F.obs=FF, p.value=pv))
    }
    if (n1<n2){
      sumy1=matrix(rep(apply(y[1:n1,],2,sum),n1),nrow=n1,byrow=T)
      sumy2=matrix(rep(apply(y,2,sum),n1),nrow=n1,byrow=T)
      dataz=x-sqrt(n1/n2)*y[1:n1,]+sqrt(1/n1*n2)*sumy1-1/n2*sumy2
      
      n=n1
      mu0=apply(x, 2, mean)-apply(y, 2, mean)
      z.bar=apply(dataz, 2, mean)
      A=(n-1)*var(dataz)
      
      T2=(n-1)*n*t(z.bar-mu0)%*%solve(A)%*%(z.bar-mu0)
      FF=(n-p)/((n-1)*p)*T2
      
      pv=round(1-pf(FF, p, n-p),3)
      return(list(zmean=z.bar, F.obs=FF, p.value=pv))
    }
    if (n1>n2){
      sumy1=matrix(rep(apply(x[1:n2,],2,sum),n2),nrow=n2,byrow=T)
      sumy2=matrix(rep(apply(x,2,sum),n2),nrow=n2,byrow=T)
      dataz=y-sqrt(n2/n1)*x[1:n2,]+sqrt(1/n2*n1)*sumy1-1/n1*sumy2
      
      n=n2
      mu0=apply(y, 2, mean)-apply(x, 2, mean)
      z.bar=apply(dataz, 2, mean)
      A=(n-1)*var(dataz)
      
      T2=(n-1)*n*t(z.bar-mu0)%*%solve(A)%*%(z.bar-mu0)
      FF=(n-p)/((n-1)*p)*T2
      
      pv=round(1-pf(FF, p, n-p),3)
      return(list(zmean=z.bar, F.obs=FF, p.value=pv))
    }
  }
}
