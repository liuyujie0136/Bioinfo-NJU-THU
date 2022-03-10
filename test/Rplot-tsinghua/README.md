# Lu lab's shared scripts for plots

output repo: https://gitlab.com/dongzhuoer/lulab-plot-output

## Overview

This book aims to provide a better way to present Lu lab's shared tips for plots. It include both bookdown and Rmarkdown version. Writen by Xiaocheng Xi, adopted by Zhuoer Dong.

## Usage

https://lulab.gitbooks.io/teaching/content/appendix/appendix.plots.html#plot-in-Rstudio

1. Install R

1. Install RStudio

1. clone the repo, double-click `lulab-plot.Rproj` in it.

1. Install the package (`Build` -> `Install and Restart`)

1. You can either open `all.Rmd`, scroll down the file and click every gree right arrow, or open `index.Rmd`, then run other `.Rmd` one by one (each file describes a plot type).

## develop

> Warnings!!!  
> Do not run the following code. It only for the maintainer.

```bash
R -e "browseURL('gitbook/index.html')"
```

```bash
cd /media/computer/work/Git/lulab-plot-output
rm -r public; cp -r ../../bookdown/lulab-plot/gitbook/ . && mv gitbook public
git add .
git cam 'update gitbook'
git push
```

## to do

- [] gitlab IDE not show preview for .Rmd


