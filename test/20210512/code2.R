library(ggplot2)
dat <- data.frame(
  year = factor(sample(2010:2014, 400, replace = T)),
  continent = factor(sample(c("EU", "US", "Asia"),
                            400, replace = T)),
  gender = factor(sample(c("male", "female"),
                         400, replace = T)),
  amount = sample(20:5000, 400, replace = T)
)

ggplot(dat, aes(
  x = paste0(continent, "_", year),
  y = amount,
  fill = paste0(continent, " ", gender)
)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  scale_x_discrete(labels = c("2010", "2011", "2012\nAsia", "2013", "2014",
                              "2010", "2011", "2012\nEU", "2013", "2014",
                              "2010", "2011", "2012\nUS", "2013", "2014")) +
  guides(fill=guide_legend(title="Gender")) +
  theme(axis.title.x = element_blank())
