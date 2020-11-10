### Draw ggplot2 graphs via clicking mouse

## find where these packages can be installed from

library(wherepackage)
d=loadData()

where(data=d,package='esquisse')
where(data=d,package='ggThemeAssist')

## install these packages

install.packages('esquisse')
install.packages('ggThemeAssist')

## use them!
library(esquisse)
library(ggThemeAssist)
library(ggplot2)

# Note: they are in addins!

## Examples

# exported code from esquisse ggplot2 builder
library(ggplot2)

p=ggplot(economics) +
 aes(x = date, y = pce, colour = unemploy) +
 geom_line(size = 1L) +
 scale_color_gradient() +
 theme_minimal()

# highlight this to use ggThemeAssist
p

# adjusted code from ggThemeAssist
p + theme(panel.background = element_rect(fill = "gray72", 
    linetype = "solid"), plot.background = element_rect(fill = "antiquewhite"))

