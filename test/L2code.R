;	# Separates two commands on the same line.
A=15; B=18
A;B
A
B
A==B
A!=B
A<B
A<=B
A>=B

X=c(1,2,4,3.3,12,15,0,19)
is.na(X)
Y=c(1,2,4,3.3,12,15,19)
Z=X&Y
Z

numeric(25) 
character(20) 
N=seq(-4,4,0.1) 
N
1:10-3
1:(10-3)#
rep("Niu",7)
Y=c(1,2,4,3.3,12,15,19)
rep(Y-2,7)
matrix(0,nrow=3,ncol=4)
matrix(1:12,nrow=3,ncol=4)
matrix(1:12,nrow=3,ncol=4,byrow=F)
x=pi
round(x,2) 
round(x,4)
log(x)	#natural log
log10(x)	#log base 10
sqrt(x)	#square root
exp(x)
sin(x);cos(x);tan(x)
asin(x);acos(x);atan(x)
length(x) 
length(Y)
min(x,Y) 

X<- seq(1,10)
for (i in x) {print("hello")}

setwd("B:/Teaching/BioSta/2019BioSta")
save.image("B:/Teaching/BioSta/2019BioSta/.RData")
load("B:/Teaching/BioSta/2019BioSta/.RData")

L2data=read.csv("L2data.csv",head=T)
L2data

fix(L2data)
L2data[1]
L2data[1,]
L2data[,1:2];
L2data[,c(1,3)];
L2data$Height
subset(L2data,Gender=="M") 
L2data[1:3] #gives first 5 rows of data

L2data[,1:3] 
L2data[L2data[,3]<169,] 
L2data[3,3]
summary(L2data) 
sum(L2data$Height)
mean(L2data$Height)
var(L2data$Height)
sd(L2data$Height)
median(L2data$Height)
cor(L2data$Height,L2data$Weight)#

hist(L2data$Height)
boxplot(L2data$Height)
boxplot(L2data[L2data$Gender=="M",]$Height,L2data[L2data$Gender=="F",]$Height)

plot(L2data$Height~L2data$Weight)
attach(L2data)
plot(Weight~Height,data=L2data)
plot(L2data$Height, L2data$Weight)
