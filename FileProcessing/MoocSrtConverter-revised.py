filename=input()
fname=filename.split('/')
t=open('/'.join(fname[:-2])+'/'+fname[-3]+'.txt','a',encoding='utf-8')
while filename!='0':
    f=open(filename,'r+',encoding='utf-8')
    s=[]
    for line in f:
        if '-->' in line:
            1==1
        elif line[0] in '0123456789':
            2==2
        else:
            line=line.replace('\n',' ')                
            s.append(line)
    t.write(fname[-2])
    t.write('\n')
    t.write(fname[-1][:-4])
    t.write('\n')
    del s[0]
    for item in s:
        t.write(item)
    t.write('\n\n\n')
    f.close()
    filename=input()
    fname=filename.split('/')
t.close()
    