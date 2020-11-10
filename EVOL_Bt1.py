p=0.001
n=1
t=0.9
h=0.1
refuge=0.5
rr=1
rs=(1-(1-h)*t)*(1-refuge)+1*refuge
ss=(1-t)*(1-refuge)+1*refuge
while p<0.1:
    R=2*pow(p,2)*rr+2*p*(1-p)*rs
    S=2*pow(1-p,2)*ss+2*p*(1-p)*rs
    p=R/(R+S)
    n+=1
print(n)
