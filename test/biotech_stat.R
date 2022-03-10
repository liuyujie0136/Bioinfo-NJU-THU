library(ggplot2)

x <- read.csv("stat.csv")

var.test(x[x$gender == 0, 4], x[x$gender == 1, 4])

t.test(x[x$gender == 0, 4], x[x$gender == 1, 4], var.equal = T)

ggplot() + geom_histogram(data = x[x$gender == 0, ],
                          mapping = aes(x = age1))

ggplot() + geom_histogram(data = x[x$gender == 0, ],
                          mapping = aes(x = age2))

ggplot() + geom_histogram(data = x[x$gender == 1, ],
                          mapping = aes(x = age1))

ggplot() + geom_histogram(data = x[x$gender == 1, ],
                          mapping = aes(x = age2))

shapiro.test(x[x$gender == 0, 4])
shapiro.test(x[x$gender == 1, 4])

wilcox.test(x[x$gender == 0, 4], x[x$gender == 1, 4])

