### Export plots to PPT

## Use "export" (https://github.com/tomwenseleers/export)

library(export)

# get system time for suitable file name
t=as.character(Sys.time())
t=gsub(" ","-",t)
t=gsub(":","-",t)
fname=paste0("Rplot-",t,".pptx")

# export plots
graph2ppt(file=fname)

# in one line (easy to copy)
export::graph2ppt(file=paste0("Rplot-",gsub(":","-",gsub(" ","-",as.character(Sys.time()))),".pptx"))


## Use "eoffice" to export graph/table to MS Office

library(ggplot2)
library(ggplotify) # convert plot to ggplot object
library(eoffice) # export

f = "eoffice.pptx"
p = as.ggplot(~plot(cars, cex.lab=2, cex.main=2,
                    xlab="biobabble", ylab="biobabble",
                    main = "Example"))

g = ggplot(data = cars) + geom_point(mapping = aes(x = speed, y = dist)) + labs(x = "biobabble", y = "biobabble", title = "Example")

topptx(p, f)
topptx(g, f)

# open file in R
library(rvcheck)
o(f)
