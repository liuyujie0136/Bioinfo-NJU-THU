import wordcloud

f = open("C:/Users/10784/Documents/Python/word-cloud.txt","r", encoding="utf-8")
txt = f.read()
f.close()

txt = txt.replace("\n"," ")

w = wordcloud.WordCloud(height=800,width=800,background_color="white",max_words=100)

w.generate(txt)

w.to_file("C:/Users/10784/Documents/Python/word-cloud.png")
