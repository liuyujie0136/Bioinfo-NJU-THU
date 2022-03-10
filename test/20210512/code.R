library(ggplot2)
library(ggsignif)

df <- read.csv("boxplot_example_1.csv")

ggplot(data = df, mapping = aes(x = paste0(dose, "_", supp), y = len)) +
  geom_violin(mapping = aes(fill = supp)) +
  geom_boxplot(width = 0.2,
               position = position_dodge(0.9)) +
  theme_bw() +
  xlab("dose") +
  scale_x_discrete(labels = c("dose1",
                              "",
                              "dose2",
                              "",
                              "dose3",
                              "")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = -1),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("yellowgreen", "violetred1")) +
  geom_signif(
    comparisons = list(c("dose_1_OJ", "dose_1_VC"),
                       c("dose_2_OJ", "dose_2_VC")),
    map_signif_level = TRUE,
    y_position = c(22, 28),
    tip_length = c(0.02, 0.36, 0.02, 0.2)
  )
