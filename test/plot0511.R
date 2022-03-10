library(ggplot2)
library(ggpubr)

WT <- c(86.6, 87.8, 93.9, 83.5, 88.6, 83.7, 88.9)
CAG_OX <- c(69.8, 92.0, 75.5, 72.8, 82.7, 71.2, 74.2, 85.3, 84.9)

data <-
  data.frame(type = c(rep("WT", length(WT)), rep("CAG-OX", length(CAG_OX))),
             value = c(WT, CAG_OX))

p <-
  ggbarplot(
    data = data,
    x = "type",
    y = "value",
    add = c("mean_sd")
  ) +
  geom_jitter(
    mapping = aes(color = type),
    size = 2,
    alpha = 0.6,
    width = 0.1,
    height = 0.1
  ) +
  tinyfuncr::geom_signif_wrapper(
    comparisons = list(c("WT", "CAG-OX")),
    test = "t.test",
    margin_top = 10,
    tip_length = 5,
    textsize = 8
  ) +
  coord_cartesian(ylim = c(0, 110)) +
  scale_y_continuous(breaks = seq(0, 200, 10), expand = c(0, 0)) +
  ylab("Preference (%)") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

ggsave("plot0511.pdf", p, height = 4, width = 4)
