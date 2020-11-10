f=open('[1]--需要：情绪和情感的催化剂.txt','r+',encoding='utf-8')
t=open('[1]--需要：情绪和情感的催化剂-revised.txt','w',encoding='utf-8')
s=[]
for line in f.readlines():
    line=line.replace('\n','')
    for ch in '0123456789:,->':
        line=line.replace(ch,'')
    s.append(line)
for item in s:
    t.write(item)
f.close()
t.close()
