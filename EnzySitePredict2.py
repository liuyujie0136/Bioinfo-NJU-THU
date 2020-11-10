f=open('EnzySitePredict_Seq.txt','r')
l53,l35,c=[],[],0
for line in f:
    line=line.replace('\n','')
    if c!=0:
        l53.append(line.upper())
    else:
        c=1
        print('Input Seq:')
        print(line)
s53=''.join(l53)
print(s53)
pair={'A':'T','T':'A','C':'G','G':'C'}
for item in s53:
    l35+=pair[item]
s35=''.join(l35)
ans={}
for i in range(len(s53)-5):
    if s53[i:i+6]==s35[i+5:i-1:-1]:
        ans[i+1]=s53[i:i+6]
if len(ans)!=0:
    print('\nPredicted Enzyme Cutting Sites:')
    for item in ans:
        print(str(item)+'#:',ans[item])
else:
    print('\nEnzyme Cutting Sites Not Found')

f.close()
print('\nPress <Enter> to quit...',end='')
input()
