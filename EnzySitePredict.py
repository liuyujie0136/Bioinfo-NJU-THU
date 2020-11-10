s53=input('Input Seq:')
s53=s53.upper()
pair={'A':'T','T':'A','C':'G','G':'C'}
l35=[]
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

print('\nPress <Enter> to quit...',end='')
input()
