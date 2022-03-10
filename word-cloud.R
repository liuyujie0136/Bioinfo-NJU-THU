## Wordcloud Visualization of TCGA Research Network Publications in CNS

library("tm")
library("SnowballC")
library("wordcloud")
library("RColorBrewer")

file <- readLines("word-cloud.txt", encoding = "UTF-8")
text <-
  data.frame(
    doc_id = 1:length(file),
    text = file,
    stringsAsFactors = FALSE
  )

# Load the data as a corpus
docs <- Corpus(DataframeSource(text))

# Replace some characters to space
toSpace <-
  content_transformer(function(x, pattern)
    gsub(pattern, " ", x))
docs <- tm_map(docs, toSpace, "/")
docs <- tm_map(docs, toSpace, "@")
docs <- tm_map(docs, toSpace, "\\|")

# Convert the text to lower case
docs <- tm_map(docs, content_transformer(tolower))

# Remove numbers
docs <- tm_map(docs, removeNumbers)

# Remove English common stopwords
docs <- tm_map(docs, removeWords, stopwords("english"))

# Remove your own stop word
#docs <- tm_map(docs, removeWords, c("a", "b"))

# Remove punctuations
docs <- tm_map(docs, removePunctuation)

# Eliminate extra white spaces
docs <- tm_map(docs, stripWhitespace)

# Text stemming
#docs <- tm_map(docs, stemDocument)

# The result
dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m), decreasing = TRUE)
d <- data.frame(word = names(v), freq = v)
head(d, 10)

# WordCloud plot
set.seed(1234)
pdf("word-cloud.pdf", height = 6, width = 6)
wordcloud(
  words = d$word,
  freq = d$freq,
  min.freq = 1,
  max.words = 200,
  random.order = FALSE,
  rot.per = 0.35,
  colors = brewer.pal(8, "Dark2")
)
dev.off()
