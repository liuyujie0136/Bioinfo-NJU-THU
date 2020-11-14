### using images and emoji in ggplot2

## package 'ggimage'

install.packages("ggimage") # OR: devtools::install_github("GuangchuangYu/ggimage")

library(ggplot2)
library(ggimage)

imgfile='C:/Users/10784/Documents/GitHub/liuyujie0136.github.io/logo.jpg'

d=data.frame(x=1:5,y=6:10,a=11:15,b=16:20,img=imgfile,size=10,replace=T)

ggplot(d)+geom_point(aes(x,y))+geom_image(aes(a,b,image=img))



## package 'emojifont'

install.packages("emojifont")

library(ggplot2)
library(emojifont)

# emoji characters

search_emoji('smile')
message(emoji(search_emoji('smile')))

# base plot

set.seed(123)
x <- rnorm(10)
set.seed(321)
y <- rnorm(10)
plot(x, y, cex=0)
text(x, y, labels=emoji('cow'), cex=1.5, col='steelblue', family='EmojiOne')

# ggplot2

d <- data.frame(x=x, y=y, label = sample(c(emoji('cow'), emoji('camel')), 10, replace=TRUE),type = sample(LETTERS[1:3], 10, replace=TRUE))
ggplot(d, aes(x, y, color=type, label=label)) +
  geom_text(family="EmojiOne", size=6)

# geom_emoji layer for easy use in ggplot2

ggplot() + geom_emoji("rose", color='steelblue') + theme_void()

x = seq(0, 2*pi, length=30)
y = sin(x)
ggplot() + geom_emoji('heartbeat', x=x, y=y, size=10)

# Font Awesome

set.seed(2016-03-09)
fa <- fontawesome(c('fa-github', 'fa-weibo', 'fa-twitter', 'fa-android', 'fa-coffee'))
d <- data.frame(x=rnorm(20), y=rnorm(20), label=sample(fa, 20, replace=T))
ggplot(d, aes(x, y, color=label, label=label)) +
  geom_text(family='fontawesome-webfont', size=6)+
  xlab(NULL)+ylab(NULL) +
  theme(legend.text=element_text(family='fontawesome-webfont'))

# geom_fontawesome layer

ggplot() + geom_fontawesome("fa-github", color='black') + theme_void()

