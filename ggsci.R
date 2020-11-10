### ggsci - Scientific Journal and Sci-Fi Themed Color Palettes for ggplot2

## short example - in Non-ggplot2 Graphics
library(ggsci)
vignette("ggsci")  # help page

color=pal_lancet("lanonc",alpha = 0.6)(4) # 4 means choose 4 colors from all
color
scales::show_col(color)

# use in ggplot2
#scale_color_<palname>()
#scale_fill_<palname>(), in this example, use scale_fill_lancet()


## Discrete Color Palettes examples

# data
library("ggsci")
library("ggplot2")
library("gridExtra")

data("diamonds")

p1 = ggplot(subset(diamonds, carat >= 2.2),
            aes(x = table, y = price, colour = cut)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", alpha = 0.05, size = 1, span = 1) +
  theme_bw()

p2 = ggplot(subset(diamonds, carat > 2.2 & depth > 55 & depth < 70),
            aes(x = depth, fill = cut)) +
  geom_histogram(colour = "black", binwidth = 1, position = "dodge") +
  theme_bw()

# Example drawing palettes
# NPG (Nature Publishing Group) palette
p1_npg = p1 + scale_color_npg()
p2_npg = p2 + scale_fill_npg()
grid.arrange(p1_npg, p2_npg, ncol = 2)

# AAAS (American Association for the Advancement of Science) palette
p1_aaas = p1 + scale_color_aaas()
p2_aaas = p2 + scale_fill_aaas()
grid.arrange(p1_aaas, p2_aaas, ncol = 2)

# UCSCGB (UCSC Genome Browser) palette
pdf('plot.pdf', width = 15, height = 8)
p1_ucscgb = p1 + scale_color_ucscgb()
p2_ucscgb = p2 + scale_fill_ucscgb()
grid.arrange(p1_ucscgb, p2_ucscgb, ncol = 2)
dev.off()


## Continuous Color Palettes examples

# data
library("reshape2")

data("mtcars")
cor = cor(unname(cbind(mtcars, mtcars, mtcars, mtcars)))
cor_melt = melt(cor)

p3 = ggplot(cor_melt,
            aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(colour = "black", size = 0.3) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

# GSEA GenePattern palette
p3_gsea = p3 + scale_fill_gsea()
p3_gsea_inv = p3 + scale_fill_gsea(reverse = TRUE)
grid.arrange(p3_gsea, p3_gsea_inv, ncol = 2)

# Material Design Color Guidelines palettes
#We generate a random data matrix first:
library("reshape2")

set.seed(42)
k = 9
x = diag(k)
x[upper.tri(x)] = runif(sum(1:(k - 1)), 0, 1)
x_melt = melt(x)

p4 = ggplot(x_melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(colour = "black", size = 0.3) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() + theme(
    legend.position = "none", plot.background = element_blank(),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_blank(), axis.text.y = element_blank(),
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    panel.background = element_blank(), panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Plot the matrix with the 19 material design color palettes:
grid.arrange(
    p4 + scale_fill_material("red"),         p4 + scale_fill_material("pink"),
    p4 + scale_fill_material("purple"),      p4 + scale_fill_material("deep-purple"),
    p4 + scale_fill_material("indigo"),      p4 + scale_fill_material("blue"),
    p4 + scale_fill_material("light-blue"),  p4 + scale_fill_material("cyan"),
    p4 + scale_fill_material("teal"),        p4 + scale_fill_material("green"),
    p4 + scale_fill_material("light-green"), p4 + scale_fill_material("lime"),
    p4 + scale_fill_material("yellow"),      p4 + scale_fill_material("amber"),
    p4 + scale_fill_material("orange"),      p4 + scale_fill_material("deep-orange"),
    p4 + scale_fill_material("brown"),       p4 + scale_fill_material("grey"),
    p4 + scale_fill_material("blue-grey"),
    ncol = 6)


### OTHER TOOLS ABOUT COLOR IN R
### Rcolorbrewer - Creates nice looking color palettes especially for thematic maps

library(RColorBrewer)

brewer.pal.info
display.brewer.pal(11, "BrBG")  # usage: display.brewer.pal(n, name)
# note:
 # 11: Number of different colors in the palette, minimum 3, maximum depending on palette
 # "BrBG": A palette name shown in brewer.pal.info
display.brewer.all()

# Make color palettes from ColorBrewer available as R palettes
color1 <- brewer.pal(8, "Set1")
color2 <- brewer.pal(8, "Dark2")

scales::show_col(color1)
scales::show_col(color2)
