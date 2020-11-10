## Export plots to PPT
# Reference: https://github.com/tomwenseleers/export

# library package "export"
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
