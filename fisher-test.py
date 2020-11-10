from scipy.special import comb
sum=0
for i in range(6,51):
    sum+=comb(208,i)*comb(21042-208,50-i)/comb(21042,50)
print(sum)

# may have bugs
