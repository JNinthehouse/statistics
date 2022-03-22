mu.test.multinorm=function(data, mu0, Sigma0=-1, alpha=0.05)   
###############################################################
## H0: mu=mu0 
## Chisq testing (Sigma is unknown)
## F testing (Sigma is known)
##############  Input  ########################################
## data  = design matrix with the ith sample in the ith line
## mu0   = mu0 for null hypothesis
## Sigma0= the known variance matrix
## alpha = the significant level, default value = 0.05
############## Output  ########################################
## Reject.area = reject region
## p.value     = p value
###############################################################
{

data=as.matrix(data)
n=nrow(data)
p=ncol(data)

if (Sigma0!=-1){
X.bar=apply(data, 2, mean)
T1=n*t(X.bar-mu0)%*%solve(Sigma0)%*%(X.bar-mu0)

a2=qchisq(1-alpha, p)

reject=matrix(c(T1, a2), nrow=1)
rownames(reject)=c("Reject")
colnames(reject)=c("Obs", ">1-alpha")

pv=1-pchisq(T1, p)
return(list(Reject.area=reject, Xmean=X.bar, Chi.obs=T1, p.value=pv))}

else if (Sigma0==-1){
X.bar=apply(data, 2, mean)
A=(n-1)*var(data)

T2=(n-1)*n*t(X.bar-mu0)%*%solve(A)%*%(X.bar-mu0)
FF=(n-p)/((n-1)*p)*T2

pv=1-pf(FF, p, n-p)
return(list(Xmean=X.bar, F.obs=FF, p.value=pv))
}
}
