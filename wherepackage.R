## find where the packages can be installed from

# remotes::install_github('yikeshu0611/wherepackage')

library(wherepackage)

d=loadData()

where(data=d,package='export')

# remotes::install_version('export','0.2.2')

# Note: run this code in RStudio, in vscode may encounter errors

# Also, https://rdrr.io/ is a good website for R Package Documentation

# Note: may have errors!