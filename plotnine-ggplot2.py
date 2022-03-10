from plotnine import ggplot, geom_jitter, ggsave, aes
from plotnine.data import mpg

mpg.head()

p=ggplot(data=mpg)+geom_jitter(aes(y='cty',x='displ',color="cyl"))

ggsave(p, "1.pdf", height=6, width=6)
