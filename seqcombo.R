### DNA sequence recombination visualization - using seqcombo

## Installation
repo = c("https://cran.rstudio.com",
         "https://bioconductor.org/packages/release/bioc",
         "https://guangchuangyu.github.io/drat")
install.packages("seqcombo", repo=repo)
#remotes::install_github('YuLab-SMU/seqcombo')


## example
library(seqcombo)
fas <- list.files(system.file("examples","GVariation", package="seqcombo"),
                  pattern="fas", full.names=TRUE)
x <- lapply(fas, seqdiff)
plts <- lapply(x, plot)
plot_grid(plotlist=plts, ncol=1, labels=LETTERS[1:3])