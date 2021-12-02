# Import Data and Required Libraries
library(ggplot2)
library(forecast)
library(tseries)
library(gridExtra)
library(MASS)
library(lmtest)


sap <- read.csv("C:/Users/Saberi/Desktop/df.csv")
sap$week<-(as.Date(sap$week,format="%m/%d/%Y"))
names(sap)<-c("week","y")
week<-sap$week[1:120]
y<-sap$y[1:120]
df<-data.frame(week,y) 


VAR.MEAN<-function(DATA)
{
  MEAN<-mean(DATA)
  n<-length(DATA)
  L<-floor(sqrt(n))
  OUT1<-acf(DATA,plot=F,type="cov")$acf[1:(L+1)]
  VARMEAN<-round(((OUT1[1]+2*(sum(OUT1[2:(L+1)]*(1-seq(1:L)/L))))/n),10)
  Z<-MEAN/sqrt(VARMEAN)
  
  return(list(MEAN=MEAN,VARMEAN=VARMEAN,Z=Z))
}
WT<-function(DATA,phi){
  wt<-0
  p<-length(phi)
  n<-length(DATA)
  for(i in (p+1):n){
    wt[i]<-DATA[i]-DATA[(i-1):(i-p)]%*%phi
  }
  Wt<-wt[(p+1):n]
  return(Wt)
}
SEACF<-function(DATA,K,L){
  
  RHOKJ<-matrix((rep(0,(K+1)*(L+1))),nrow=K+1)
  TESTRHOKJ<-matrix((rep(0,(K+1)*(L+1))),nrow=K+1)
  TEST<-matrix((rep(0,(K+1)*(L+1))),nrow=K+1)
  n<-length(DATA)
  DATA<-DATA-mean(DATA)
  phihat<-0
  Data<-matrix(rep(0,n*(K+1)),ncol=(K+1))
  
  for(i in 1:(K+1)){Data[1:(n-i+1),i]<-DATA[i:n]}
  Err<-matrix(rep(0,n*(L+2)),ncol=L+2)
  
  for(k in 1:K){				
    Out<-lm(Data[1:(n-k),(k+1)]~Data[1:(n-k),1:k])
    
    Err[1:n,1]<-c(Out$resid,rep(0,k))
    
    Out<-lm(Data[2:(n-k),(k+1)]~Data[2:(n-k),1:k]+Err[1:(n-k-1),1])
    phihat[1:k]<-Out$coeff[2:(k+1)]
    
    RHOKJ[k+1,1]<-acf(WT(DATA,phihat[k:1]),plot=F)$acf[2]
    Err[1:n,2]<-c(Out$resid,rep(0,k+1))
    
    for(l in 2:(L+1)){
      
      Out<-lm(Data[(l+1):(n-k),(k+1)]~Data[(l+1):(n-k),1:k]+Err[1:(n-k-l),l:1])
      
      Err[1:n,(l+1)]<-c(Out$resid,rep(0,k+l))
      phihat[1:k]<-Out$coeff[2:(k+1)]
      RHOKJ[k+1,l]<-acf(WT(DATA,phihat[k:1]),plot=F,lag=l+1)$acf[l+1]		
    }
    Err<-matrix(rep(0,n*(L+2)),ncol=L+2) }
  RHOKJ[1,]<-acf(DATA,plot=F,lag=L+2)$acf[2:(L+2)]
  RHOKJ<-round(RHOKJ,2)
  
  for(i in 1:(K+1)){
    for(j in 1:(L+1)){
      if(abs(RHOKJ[i,j])> 2*sqrt((n-i-j)^(-1)))
        TESTRHOKJ[i,j]<-1
      TEST[i,j]<-round(abs(RHOKJ[i,j])/(2*sqrt((n-i-j)^(-1))),2)
    }
  }
  dimnames(RHOKJ)<-list(paste("p=",0:K),paste("q=",0:L))
  dimnames(TESTRHOKJ)<-list(paste("p=",0:K),paste("q=",0:L))
  
  return(list(RHOKJ=RHOKJ,TESTRHOKJ=TESTRHOKJ,TEST=TEST))
}

PORT.TEST<-function(N, OUT, K)
{
  tilRHO1 <- matrix(rep(0, (K + 1) * (K + 1)),	ncol = K + 1)
  
  acfer <- acf(OUT,plot=F,type="cor")$acf[1:(K + 1)]
  pacfer <- acf(OUT,plot=F,type="pa")$acf[1:K]
  
  
  QBP <- sum(acfer[2:(K + 1)]^2) * N
  QLB <- N * (N + 2) * sum(acfer[2:(K + 1)]^2 *(N - 1:K)^(-1))
  QMT <- N * (N + 2) * sum(pacfer[1:K]^2 * (N -1:K)^(-1))
  
  for(i in 1:(K + 1)) {
    for(j in 1:(K + 1)) {
      tilRHO1[i, j] <- acfer[abs(i-j) + 1]
    }
  }
  
  HATD <- N * (1 - (prod(eigen(tilRHO1)$values	))^(1/K))
  
  return(list(QBP=QBP, QLB=QLB, QMT=QMT, HATD=HATD))
}


# Be Friend With Data
theme_set(theme_bw())
ggplot(data = df, aes(x = week, y = y))+
  geom_line(color="#00AFBB" , size=1)

# check wether we need variance stablzer transformations or not
x<-(1:length(df$y))
data<-df$y
m=lm(data ~ x)
boxcox(m)

# show the nonstationarity of the data by dickey fuller test
adf.test(df$y)

# making data stationary by differentiating
ndiffs(df$y,alpha = 0.05,test = c("kpss", "adf", "pp")
       ,type = c("level", "trend"),max.d = 5)

df$temp<-c(0,diff(df$y))
adf.test(df$temp)

ggplot(data = df, aes(x = week, y = temp))+
  geom_line(color = "#00AFBB",size=1)

ggtsdisplay(df$y)

ggtsdisplay(df$temp,smooth = TRUE)

# DETECT SOME SUTAIBLE CANDIDATES TO MODEL OUR DATA PROPERLY
# find a suitable candidate from moving-average models
# build thetahats for hypothesis test 
thetahat=acf(df$temp,10,type="correlation",plot = F)$acf[2:11]
thetahat

zscore<-0
zscore[1]<-sqrt(length(df$y))*thetahat[1]/sqrt(1+2*(thetahat[1])^2)

for (i in 2:10){
  zscore[i]=sqrt(length(df$y))*thetahat[i]/sqrt((1+2*sum(thetahat[1:(i-1)]^2)))
  if(abs(zscore[i])> 1.96){
    print(i)
  }
}

# find a suitable candidate from  auto regrresive models
# build phihat for hypothesis test
phihat=acf(df$temp,10,type="partial",plot = F)$acf

zzscore<-0
for (i in 1:10){
  zzscore[i]=sqrt(length(df$y))*phihat[i]
  if(abs(zzscore[i])>1.96){
    print(i)
  }
}


# find a suitable candidate from  ARMA models
SEACF(df$temp,6,6)$TESTRHOKJ

# check wether we have drift or not
VAR.MEAN(df$temp)

# detect our final models for investigating
Arima(df$y, order = c(4,1,0),include.drift = F) #from auto reggresive model
#we detect all the below models from Extended Autocorrelation Function
Arima(df$y, order = c(0,1,2),include.drift = F)

Arima(df$y, order = c(0,1,1),include.drift = F)

Arima(df$y, order = c(1,1,2),include.drift = F)

Arima(df$y, order = c(2,1,3),include.drift = F)

Arima(df$y, order = c(3,1,4),include.drift = F)

Arima(df$y, order = c(4,1,4),include.drift = F)


# model validation
fit1<-Arima(df$y, order = c(4,1,0),include.drift = F)
fit2<-Arima(df$y, order = c(0,1,2),include.drift = F)
fit3<-Arima(df$y, order = c(0,1,1),include.drift = F)
fit4<-Arima(df$y, order = c(1,1,2),include.drift = F)
fit5<-Arima(df$y, order = c(3,1,4),include.drift = F)

# investigating model 1 ARIMA(4,1,0)
checkresiduals(fit1)
tsdiag(fit1)

SEACF(fit1$residuals,6,6)$TESTRHOKJ

ks.test(fit1$residuals,mean(fit1$residuals),sd(fit1$residuals))

qqnorm(fit1$residuals)

PORT.TEST(120,fit1$residuals,15)
qchisq(0.95,14)

coeftest(fit1)

# investigating model 2 ARIMA(0,1,2)
checkresiduals(fit2)
tsdiag(fit2)

SEACF(fit2$residuals,6,6)$TESTRHOKJ

ks.test(fit2$residuals,mean(fit2$residuals),sd(fit2$residuals))

qqnorm(fit2$residuals)

PORT.TEST(120,fit2$residuals,15)
qchisq(0.95,14)

coeftest(fit2)

# investigating model 3 ARIMA(0,1,1)
checkresiduals(fit3)
tsdiag(fit3)

SEACF(fit3$residuals,6,6)$TESTRHOKJ

ks.test(fit3$residuals,mean(fit3$residuals),sd(fit3$residuals))

qqnorm(fit3$residuals)

PORT.TEST(120,fit3$residuals,15)
qchisq(0.95,14)

coeftest(fit3)

# investigating model 4 ARIMA(1,1,2)
checkresiduals(fit4)
tsdiag(fit4)

SEACF(fit4$residuals,6,6)$TESTRHOKJ

ks.test(fit4$residuals,mean(fit4$residuals),sd(fit4$residuals))

qqnorm(fit4$residuals)

PORT.TEST(120,fit4$residuals,15)
qchisq(0.95,14)

coeftest(fit4)

# investigating model 5 ARIMA(3,1,4)
checkresiduals(fit5)
tsdiag(fit5)

SEACF(fit5$residuals,6,6)$TESTRHOKJ

ks.test(fit5$residuals,mean(fit5$residuals),sd(fit5$residuals))

qqnorm(fit5$residuals)

PORT.TEST(120,fit5$residuals,15)
qchisq(0.95,14)

coeftest(fit5)

### forecast
df$fit4<-fit4$fitted
ggplot()+
  geom_line(data = df, aes(x = week, y = fit4),color = "black",size=1)+
  geom_line(data = df, aes(x = week, y = y),color = "#00AFBB",size=1)

future4 = forecast(fit4, h = 12)
plot(future4)
df$fit4<-fit4$fitted

est4<-as.data.frame(forecast(fit4,h=12,level=c(95)))
yhat4<-c(rep(NA,120),est4$`Point Forecast`)
sap<-cbind(sap,yhat4)
ggplot(sap)+
  geom_line(aes(x=week,y=y),col="royalblue3",lwd=1)+
  geom_line(aes(x=week,y=yhat4),col="deeppink3",lwd=1.3)

xhat=est4$`Point Forecast`

x=sap$y[121:132]

MSE=sum((x-xhat)*(x-xhat))/length(x)
MSE
