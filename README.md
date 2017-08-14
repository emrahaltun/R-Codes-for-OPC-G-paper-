# R-Codes-for-OPC-G-paper-
"The Odd Power Cauchy family of distributions: Properties, Regression Models and Applications"
\textbf{Developed R Code of Fitted proposed OPC-L distribution of first and second application with profile likelihood functions:}

###Profile likelihood
#pdf of OPC-Lindley distribution
library(VGAM)
OPCLd<-function(x,alpha,theta){
 G<-rep(0,0)
 for(i in 1:length(x)){
    G[i]<-(2*alpha*dlind(x[i], theta)*(plind(x[i], theta)*
    (1-plind(x[i], theta)))^(alpha-1))/(pi*((plind(x[i], theta)^(2*alpha))
    +((1-plind(x[i], theta))^(2*alpha))))
 }
return(G)
}
###App 1
library("AdequacyModel")
x<-c(0.654, 0.613, 0.315, 0.449, 0.297, 0.402, 0.379, 0.423, 0.379, 0.3235
,0.269, 0.740, 0.418, 0.412, 0.494, 0.416, 0.338, 0.392, 0.484, 0.265)
#####
cdf_OPCL<-function(par,x){
  alpha=par[1]
  theta=par[2]
  CDF=(2/pi)*atan((plind(x, theta)/(1-plind(x, theta)))^alpha)
  PDF=(2*alpha*dlind(x, theta)*(plind(x, theta)^(alpha-1))
  *((1-plind(x, theta))^(alpha-1)))/(pi*(plind(x, theta)^(2*alpha))
  +((1-plind(x, theta))^(2*alpha)))
return(CDF)
  }

pdf_OPCL<-function(par,x){
  alpha=par[1]
  theta=par[2]
  CDF=(2/pi)*atan((plind(x, theta)/(1-plind(x, theta)))^alpha)
  PDF=(2*alpha*dlind(x, theta)*(plind(x, theta)
  *(1-plind(x, theta)))^(alpha-1))/(pi*((plind(x, theta)^(2*alpha))
  +((1-plind(x, theta))^(2*alpha))))
return(PDF)
}
OPCL<-goodness.fit(pdf=pdf_OPCL, cdf=cdf_OPCL,starts = c(3,.9),
data = x,method="N", domain=c(0.1,Inf),mle=NULL)
####################################
########loglikelihood function
####
par(mfrow=c(1,2))
OPCLlikelihood<-function(par){
  alpha=par
  theta=OPCL$mle[2]
    G<- -sum(log(OPCLd(x,alpha,theta)))
return(G)
}
H<-rep(0,0)
Alpha<-seq(.01,30,0.1)
for(i in 1:length(Alpha)){
H[i]<--OPCLlikelihood(Alpha[i])
}
plot(Alpha,H,type="l",col="red",
ylab="profile Loglikelihood",xlab=expression(paste(alpha)))
########loglikelihood function
OPCLlikelihood<-function(par){
  alpha=OPCL$mle[1]
  theta=par
    G<- -sum(log(OPCLd(x,alpha,theta)))
return(G)
}
H<-rep(0,0)
teta<-seq(.01,30,0.1)
for(i in 1:length(teta)){
H[i]<--OPCLlikelihood(teta[i])
}
plot(teta,H,type="l",col="red",
ylab="profile Loglikelihood",xlab=expression(paste(theta)))
#############################################
#############################################
###App 2
x<-c(1.1,1.4,1.3,1.7,1.9,1.8,1.6,2.2,1.7,2.7,
4.1,1.8,1.5,1.2,1.4,3,1.7,2.3,1.6,2)
OPCL<-goodness.fit(pdf=pdf_OPCL, cdf=cdf_OPCL,starts = c(3,.8),
data = x,method="N", domain=c(0.1,Inf),mle=NULL)
####################################
########loglikelihood function
par(mfrow=c(1,2))
OPCLlikelihood<-function(par){
  alpha=par
  theta=OPCL$mle[2]
    G<- -sum(log(OPCLd(x,alpha,theta)))
return(G)
}
H<-rep(0,0)
Alpha<-seq(.01,30,0.1)
for(i in 1:length(Alpha)){
H[i]<--OPCLlikelihood(Alpha[i])
}
plot(Alpha,H,type="l",col="red",
ylab="profile Loglikelihood",xlab=expression(paste(alpha)))
########loglikelihood function
OPCLlikelihood<-function(par){
  alpha=OPCL$mle[1]
  theta=par
    G<- -sum(log(OPCLd(x,alpha,theta)))
return(G)
}
H<-rep(0,0)
teta<-seq(.01,10,0.1)
for(i in 1:length(teta)){
H[i]<--OPCLlikelihood(teta[i])
}
plot(teta,H,type="l",col="red",
ylab="profile Loglikelihood",xlab=expression(paste(theta)))


\textbf{Developed R Code for Graphical Simulation}

library(AdequacyModel)

cdfnolln=function(par,x)
{
a=par[1]
mu=par[2]
sigma=par[3]
G=pnorm(x,mu,sigma)
g=dnorm(x,mu,sigma)

f=(2/pi)*atan((G/(1-G))^a)

return(f)
}

pdfnolln=function(par,x)
{
a=par[1]
mu=par[2]
sigma=par[3]
G=pnorm(x,mu,sigma)
g=dnorm(x,mu,sigma)

f=2*a*g*G^(a-1)*(1-G)^(a-1)/(pi*(G^(2*a)+(1-G)^(2*a)))

return(f)
}

qq=function(a,mu,sigma,u) {

qnorm((tan(pi*u/2))^(1/a)/(1+(tan(pi*u/2))^(1/a)),mu,sigma)

}

for( k in 1:191)
{
for( j in 1:1000)
{
datt=NULL

for(i in 1:(45+k*5)) {

p=runif(1,0,1)

a1=0.5
mu1=0.5
sigma1=2


datt[i]=qq(a1,mu1,sigma1,p)
}

fnoll=goodness.fit(pdf=pdfnolln, cdf=cdfnolln,
starts = c(a1,mu1,sigma1), data = datt,
method="N", domain=c(0,Inf))

lambda[j]=fnoll$mle[1]
beta[j]=fnoll$mle[2]
a[j]=fnoll$mle[3]

slambda[j]=fnoll$Erro[1]
sbeta[j]=fnoll$Erro[2]
sa[j]=fnoll$Erro[3]

llambda[j]=lambda[j]-1.96*slambda[j]
ulambda[j]=lambda[j]+1.96*slambda[j]
llbeta[j]=beta[j]-1.96*sbeta[j]
ubeta[j]=beta[j]+1.96*sbeta[j]
la[j]=a[j]-1.96*sa[j]
ua[j]=a[j]+1.96*sa[j]


}

ALlambda[k]=(3.92/NROW((slambda)))*sum((slambda))
ALbeta[k]=(3.92/NROW((sbeta)))*sum((sbeta))
ALa[k]=(3.92/NROW((sa)))*sum((sa))

cplambda[k]=mean((llambda) < a1 & (ulambda) > a1);
cpbeta[k]=mean((llbeta) < mu1 & (ubeta) > mu1);
cpa[k]=mean((la) < sigma1 & (ua) >sigma1);

biaslambda[k]=mean(lambda-a1)
biasbeta[k]=mean(beta-mu1)
biasa[k]=mean(a-sigma1)


mselambda[k]=sum((lambda-a1)^2)/1000
msebeta[k]=sum((beta-mu1)^2)/1000
msea[k]=sum((a-sigma1)^2)/1000

}


\textbf{Developed R Code for LOPC-W regression model}

sy=function(par)
{
alpha=par[1]
sigma=par[2]
beta0=par[3]
beta1=par[4]

a=(1/sigma)
b=exp(beta0*xxx+beta1*xx)

G=1-exp(-(y/b)^a)
g=((a/b)*(y/b)^(a-1)*exp(-(y/b)^a))


f=1-(2/pi)*atan((G/(1-G))^alpha)
return(f)
}

fx=function(par)
{

alpha=par[1]
sigma=par[2]
beta0=par[3]
beta1=par[4]

a=(1/sigma)
b=exp(beta0*xxx+beta1*xx)
G=1-exp(-(y/b)^a)
g=((a/b)*(y/b)^(a-1)*exp(-(y/b)^a))


f=2*alpha*g*G^(alpha-1)*(1-G)^(alpha-1)/(pi*(G^(2*alpha)+(1-G)^(2*alpha)))
return(f)
}

logp=function(par) {

-sum((censur)*log(fx(par)))-sum((1-censur)*log(sy(par)))

}
init=c(1,1.071,3.003,-1.052)
f=optim(par=init,fn=logp,method="BFGS",hessian=T)
