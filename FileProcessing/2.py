filename=input()
filename=filename.replace('"','')
while filename!='0':
    f=open(filename,'r+',encoding='utf-8')
    t=open(filename[:-4]+'.txt','w',encoding='utf-8')
    s=[]
    for line in f:
        line=line.replace('\n','')
        for ch in '0123456789:,->':
            line=line.replace(ch,'')
        s.append(line)
    for item in s:
        t.write(item)
    f.close()
    t.close()
    filename=input()

    
