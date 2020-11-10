### Get You to Your Desired Chart Plot Faster - using ggcharts

install.packages("ggcharts")
library(ggcharts)

## An Quick Example
# Using ggplot2 - a lot of code!
library(dplyr)
library(ggplot2)
library(ggcharts)
data("biomedicalrevenue")
biomedicalrevenue %>%
  filter(year %in% c(2012, 2015, 2018)) %>%
  group_by(year) %>%
  top_n(10, revenue) %>%
  ungroup() %>%
  mutate(company = tidytext::reorder_within(company, revenue, year)) %>%
  ggplot(aes(company, revenue)) +
  geom_col() +
  coord_flip() +
  tidytext::scale_x_reordered() +
  facet_wrap(vars(year), scales = "free_y")

# Bar chart using ggcharts - tidy and beautiful!
biomedicalrevenue %>%
  filter(year %in% c(2012, 2015, 2018)) %>%
  bar_chart(x = company, y = revenue, facet = year, top_n = 10)


## Charts Gallery

# biomedicalrevenue
data("biomedicalrevenue")
biomedicalrevenue %>%
  filter(year == 2018) %>%
  lollipop_chart(x = company, y = revenue, threshold = 30) +
  labs(
    x = NULL,
    y = "Revenue",
    title = "Biomedical Companies with Revenue > $30Bn."
  ) +
  scale_y_continuous(
    labels = function(x) paste0("$", x, "Bn."),
    expand = expansion(mult = c(0, .05))
  )

# popeurope
data("popeurope")
dumbbell_chart(
  data = popeurope,
  x = country,
  y1 = pop1952,
  y2 = pop2007,
  top_n = 10,
  point_colors = c("lightgray", "#494F5C")
) +
  labs(
    x = NULL,
    y = "Population",
    title = "Europe's Largest Countries by Population in 2007"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    labels = function(x) paste(x, "Mn.")
  )

# mtcars
data(mtcars)
mtcars_z <- dplyr::transmute(
  .data = mtcars,
  model = row.names(mtcars),
  hpz = scale(hp)
)

diverging_bar_chart(data = mtcars_z, x = model, y = hpz)

diverging_lollipop_chart(
  data = mtcars_z,
  x = model,
  y = hpz,
  lollipop_colors = c("#006400", "#b32134"),
  text_color = c("#006400", "#b32134")
)

# popch
data("popch")
pyramid_chart(data = popch, x = age, y = pop, group = sex)


## Themes Gallery
ggcharts_set_theme("theme_hermit")
bar_chart(data = diamonds, x = cut)

ggcharts_set_theme("theme_ng")
bar_chart(data = diamonds, x = cut)

ggcharts_set_theme("theme_nightblue")
bar_chart(data = diamonds, x = cut)

