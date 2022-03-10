
save.image("C:/Users/10784/Documents/R/.RData")
load("C:/Users/10784/Documents/R/.RData")


gender=c("male","female","male","female","male","female","male","female","male","female","male","female","male","female")
age=c(25,26,26,27,27,27,28,28,29,30,45,45,56,57)
weight=c(120,124,114,134,126,134,178,167,160,110,220,165,156,165)

Wdata=data.frame(gender,age,weight)
Wdata
fix(Wdata)

write.table(Wdata, file="Wdata.txt", sep="\t")

Wdata=read.table("Wdata.txt",h=T,sep="\t");
fix(Wdata)

par('mar'=c(2,2,2,2))
barplot(Wdata$weight,col="red")

library(ggplot2)

ggplot(Wdata, aes(x=gender,y=weight, fill = gender,col=gender)) + 
  geom_bar(stat="identity", position = "dodge") 


Y= c(18, 12,4, 12, 8)
X =c("US", "UK", "Australia", "Germany", "France")
pie(Y, labels = X, main="Pie Chart of Countries")

Y1=round(Y/sum(Y)*100)
X1=paste(X, Y1) # add percents to labels
X2 =paste(X1,"%",sep="") # ad % to labels
pie(Y, labels = X2, col=rainbow(length(X2)),main="Countries Pie ")

library(plotrix)
pie3D(Y, labels=X2,explode=0.1,main="Countries Pie")


Score=c(99,80,85,91,92,95,95,96,92,78,95,87,83,82,76,67,81,90,93,78,83,73,74,87,88,90,89,76,84,98,86,87,85,82)
par('mar'=c(2,2,2,2))
boxplot(Score)
mean(Score)

### For two or more group 
Wdata
aggregate(weight~gender, data=Wdata,FUN=mean)
tapply(Wdata$weight,Wdata$gender, mean)

mode(Score)
as.numeric(names(which.max(table(Score))))

statmod =function(x) {
  z =table(as.vector(x))
  names(z)[z == max(z)]
}
statmod(Score)


library(ggplot2)
ggplot(data=Wdata,aes(x=gender,y=weight,col=gender))+
  geom_boxplot(aes(fill=weight))+
  geom_point(size=5) #geom_jitter()


library(sciplot);
bargraph.CI(x.factor=gender,
            response=weight,ylim=c(0,200),
            err.width=0.02,data =Wdata,
            cex.names = 1.5, col=c("green","red"))

range(Score)
var(Score)
sd(Score)
quantile(Score)
fivenum(Score)
summary(Score) 

boxplot(Score)

aggregate(weight~gender, data=Wdata,FUN=sd)
aggregate(weight~gender, data=Wdata,FUN=quantile)

tapply(Wdata$weight,Wdata$gender, sd)
tapply(Wdata$weight,Wdata$gender, quantile)
tapply(Wdata$weight,Wdata$gender, median)

boxplot(weight~gender,data=Wdata,ylab="weight of lkkkk", col=c("green","red"))

library(psych)
describeBy(Wdata,group=gender)

library(plyr)
mWdata= ddply(Wdata, c("gender"), summarise,
              N = length(weight),
              mean = mean(weight),
              sd   = sd(weight),
              se   = sd / sqrt(N)
)
mWdata
ggplot(mWdata,aes(x=gender,y=mean, col=gender))+
  geom_point(size=8)+ 
  geom_errorbar(width=0.1, aes(ymin=mean-sd, ymax=mean+sd))

ggplot(data=Wdata,
       aes(x=gender,y=weight,col=gender))+
  geom_boxplot(aes(fill=weight))+
 geom_jitter() 

skewness(Score)
kurtosis(Score)

all.moments(Score, order.max=4 )

mystats<-function(x,na.omit=F){
  if(na.omit)
    x<-x[!is.na(x)]
  m<-mean(x)
  n<-length(x)
  s<-sd(x)
  skew<-sum((x-m)^3/s^3)/n
  kurt<-sum((x-m)^4/s^4)/n-3 
  return(c(n=n,mean=m,stdev=s,skew=skew,kurtosis=kurt))
}
mystats(Score) 

describe(Score) 



