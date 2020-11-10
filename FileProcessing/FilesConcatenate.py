filename=input()
fname=filename.split('/')
t=open('/'.join(fname[:-1])+'/合并.txt','a',encoding='utf-8')
while filename!='0':
    f=open(filename,'r+',encoding='utf-8')
    s=[]
    for line in f:
        s.append(line)
    t.write(fname[-1])
    t.write('\n')
    for item in s:
        t.write(item)
    t.write('\n\n\n')
    f.close()
    filename=input()
    fname=filename.split('/')
t.close()


