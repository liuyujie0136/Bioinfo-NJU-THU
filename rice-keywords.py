import wordcloud

f = open("C:/Users/10784/Documents/Python/rice-keywords.txt",
         "r",
         encoding="utf-8")
txt = f.read()
f.close()

txt = txt.replace("\n", " ")

w = wordcloud.WordCloud(height=800,
                        width=800,
                        background_color="white",
                        max_words=60)

w.generate(txt)

w.to_file("C:/Users/10784/Documents/Python/rice-keywords.png")
