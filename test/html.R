library(tidyverse)

file <- dir(pattern = "html$")

html <- read_table(file[2], col_names = FALSE)

html <- html[-(1:5), 1]

for (i in 1:nrow(html)){
  html[[1]][i] <- str_replace(html[[1]][i], 
                              '.+HREF="(.+)" (LAST|ADD).+>(.+)</A>', 
                              "* [\\3](\\1)")
  html[[1]][i] <- str_replace(html[[1]][i], '<TITLE>(.+)</TITLE>', "# \\1")
  html[[1]][i] <- str_replace(html[[1]][i], '.+>(.+)</H3>', "## \\1")
  html[[1]][i] <- str_replace(html[[1]][i], '&amp;', "&")
  html[[1]][i] <- str_replace(html[[1]][i], 
                              '(</?DL><p>)|(<H1>Bookmarks</H1>)', "")
  html[[1]][i] <- str_replace(html[[1]][i], 
                              '.+HREF="(.+zhihu.+)" (LAST|ADD).+>(.+)<?.+',
                              "* [\\3](\\1)")
}

write.table(html, file="bookmarks.md", col.names = FALSE, row.names = FALSE, quote = FALSE)
