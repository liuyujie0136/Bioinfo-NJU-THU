### using images and emoji in ggplot2

## package 'ggimage'

install.packages("ggimage"); devtools::install_github("GuangchuangYu/ggimage")

library(ggplot2)
library(ggimage)

imgfile='C:/Users/10784/Documents/GitHub/liuyujie0136.github.io/logo.jpg'

d=data.frame(x=1:5,y=6:10,a=11:15,b=16:20,img=imgfile,size=10,replace=T)

ggplot(d)+geom_point(aes(x,y))+geom_image(aes(a,b,image=img))



## package 'emojifont'

# note: may have errors!

install.packages("emojifont"); devtools::install_github("richfitz/remoji")

library(ggplot2)
library(emojifont)
library(remoji)

# base plot

set.seed(123)
x <- rnorm(10)
set.seed(321)
y <- rnorm(10)
plot(x, y)
text(x,y,labels=emoji('cow'), cex=1.5, col='steelblue')

# ggplot2
dd=data.frame(x=emoji(c("satisfied", "disapointed")), y=c(50, 10))
emoji_text=element_text(family="OpenSansEmoji", size=20)
ggplot(dd, aes(x, y)) + geom_bar(stat='identity', aes(fill=x)) +
  ggtitle(paste(emoji(c("+1", "-1")), collapse=" "))+
  theme(axis.text.x = emoji_text, legend.text=emoji_text, title=emoji_text) +
  xlab(NULL)+ylab(NULL)
