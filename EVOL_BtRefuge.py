q1=eval(input('Resistance allelic frequency(q):'))
t=eval(input('Relative toxicity of transgenic insecticidal cultivars(TIC):'))
h=eval(input('Dominance of resistance trait:'))
for r in range(0,100):
    refuge=r/100
    rr=1.0
    rs=(1-(1-h)*t)*(1-refuge)+1*refuge
    ss=(1-t)*(1-refuge)+1*refuge
    n=1
    q=q1
    while q<0.1:
        R=2*pow(q,2)*rr+2*q*(1-q)*rs
        S=2*pow(1-q,2)*ss+2*q*(1-q)*rs
        q=R/(R+S)
        n+=1
    if n>=100:
        print('\nRecommended refuge percentage is {}% '
              'accoding to reaching q = 0.1 after 100 generations.'.format(r))
        break
print('\nAnalysis end. Press <Enter> to quit...',end='')
input()
              
