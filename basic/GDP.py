import matplotlib.pyplot as plt
f=open('GDP.csv','r')
s,x,y=[],[],[]
for line in f:
    line=line.replace('\n','')
    line=line.split(',')
    del(line[2])
    s.append(line)
xname,yname=s[0]
del(s[0])
s=s[::-1]

for item in s:
    x.append(eval(item[0]))
    y.append(eval(item[1]))

plt.plot(x,y)
plt.xlabel(xname)
plt.ylabel(yname)
plt.show()
f.close()
