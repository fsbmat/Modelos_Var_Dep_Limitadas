rm(list=ls())
cat("\014")
a<-seq(-3,4,length.out = 20000)
#a<-ifelse(b>0,b,0)
f<-function(x){
  return((1/sqrt(2*pi))*exp(-(x^2)/2))
}
integrate(f,-Inf,0)
F1<-function(x){
  return(pnorm(x))
}
IMR<-function(x){
  f(x)/(1-F1(x))}
F1(0)
IMR(0)
IMR(8.1)
plot(IMR(a),type = "l",ylim=c(-1,5))
abline(h=0)
abline(h=4.3)

d<-function(x){
  IMR(x)*(IMR(x)-x)}

plot(d(a),type = "l",ylim = c(-0.1,1.1))
abline(h=0)
abline(h=1)
mill <- function(x) {
  return((1 - pnorm(x)) / dnorm(x))
}

plot(mill(a),type = "l")


####Playing with Truncated Normals####

set.seed(42)
x<-rnorm(200,0,4)
y<-x[x>-2]

plot(density(y),xlim=c(-16,16),main="Truncated Normal Distribution",xlab="x")
lines(density(x),col=2)
segments(x0=-2,y0=0,y1=0.15,lty=2,lwd=2)
legend(6,0.125,legend=c("Original","Truncated","Truncation\nPoint"),col=c(2,1,1),lty=c(1,1,2));

#Plot CDF because we will need the original CDF for the Inverse Mills Ratio
plot(ecdf(x),col=2,main="Empirical CDF")
lines(ecdf(y))
legend(-10,0.8,legend=c("Original","Truncated"),col=c(2,1),lty=c(1,1))

#Create 1000 Random normal variates and sample repeatedly
x.source<-rnorm(1000,0,4)
s.sample<-matrix(0,200,100)
for(i in 1:100) {s.sample[,i]<-sample(x.source,200,replace=TRUE)}

#Trunc Sum is our row means for various outputs from the bootstrap
trunc.sum<-matrix(0,200,4)
trunc.sum[,1]<-rowMeans(s.sample)
for(i in 1:200) {trunc.sum[i,2]<-mean(s.sample[i,s.sample[i,]>-2])}

#Alpha is basically a scale factor.  It gives us a scale free measure of the
#truncation point
alpha<-matrix(0,200,1)
for (i in 1:200) {alpha[i]<-(-2-trunc.sum[i,1])/sd(s.sample[i,])}

#This is the whole reason I created the file.  I didn't know of the pdf/cdf used
#for the Inverse Mills ratio were standardized or not.  The original sd was 4,
#standardized would be 1
lambda4<-dnorm(alpha,0,4)/(1-pnorm(alpha,0,4))
lambda1<-dnorm(alpha,0,1)/(1-pnorm(alpha,0,1))
#Expected value of truncated distribution for each
for (i in 1:200) {trunc.sum[i,3]<-trunc.sum[i,1]+sd(s.sample[i,])*lambda4[i]}
for (i in 1:200) {trunc.sum[i,4]<-trunc.sum[i,1]+sd(s.sample[i,])*lambda1[i]}
colnames(trunc.sum)<-c("Orig.mean","trunc.mean","sd4.mean","sd1.mean")

#did this mostly because it looks purdy
plot(density(s.sample[1,]),xlim=c(-20,20),type="n",main="100 Random Samples",xlab="x",ylim=c(0,0.225))
for (i in 2:100) {lines(density(s.sample[i,]),lwd=0.2)}
for (i in 1:100) {lines(density(s.sample[i,s.sample[i,]>-2]),lwd=0.2,col=2)};
#Here is the meat of it.  I didn't do variance but because delta follows 
#mechanically from lambda I didn't need to.
plot(density(trunc.sum[,4]),type="l",xlim=c(-1.2,3.5),main="Comparison of Bootstraped Truncated Means",xlab="Conditional Mean")
lines(density(trunc.sum[,2]),col=3)
lines(density(trunc.sum[,3]),col=4)
trunc.leg<-c("True Truncated Mean","Computed Mean\n with Standardization\n","Computed Mean\n without Standardization")
legend(x=-1.25,y=1.39,legend=trunc.leg,col=c(1,3,4),lty=c(1,1,1),cex=0.65,bty="n");
