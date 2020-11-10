import random
import re
bases = ['A','T','C','G']
SequenceList = []
Length_of_Seq=2500000
for n in range(Length_of_Seq):
    SequenceList.append(random.choice(bases))
Seq = ''.join(SequenceList)
SearchPattern1 = 'TATAAT'
SearchPattern2 = 'TTGA.A...................TA...T'
SearchPattern3 = 'TTGACA.{15,25}TATAAT'

result1 = re.finditer(SearchPattern1,Seq)
result2 = re.finditer(SearchPattern2,Seq)
result3 = re.finditer(SearchPattern3,Seq)

nmatches = 0
for match in result1:
    nmatches += 1
    color='green'
    print('Search Pattern1 finds:',match.group(),'Start Position:',match.start(),'End Position:',match.end())
print('Number of search hits of SearchPattern1= ',nmatches)
nmatches = 0
for match in result2:
    nmatches += 1
    print('Search Pattern2 finds:',match.group(),'Start Position:',match.start(),'End Position:',match.end())
print('Number of search hits of SearchPattern2= ',nmatches)
nmatches = 0
for match in result3:
    nmatches += 1
    print('Search Pattern3 finds:',match.group(),'Start Position:',match.start(),'End Position:',match.end())
print('Number of search hits of SearchPattern3= ',nmatches)

print('\nPress <Enter> to quit...',end='')
input()
