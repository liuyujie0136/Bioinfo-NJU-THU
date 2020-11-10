import jieba
f=open('沉默的羔羊.txt','r',encoding='utf-8')
txt=f.read()
words=jieba.lcut(txt)
count={}
for word in words:
    if len(word)==1:
        continue
    else:
        count[word]=count.get(word,0)+1
ls=list(count.items())
ls.sort(key=lambda x:x[1], reverse=True)
m=ls[0][0]
print(m)





