import time
import math
scale = 50
x=0
print("执行开始".center(scale//2, "-"))
start = time.perf_counter()
for i in range(scale+1):
    a = '#' * i
    b = '.' * (scale - i)
    c = (i/scale)*100
    dur = time.perf_counter() - start
    print("\r{:^3.0f}%[{}{}]{:.2f}s".format(c,a,b,dur),end='')
    time.sleep(x+(1-math.sin(x*3.14*2+3.14/2)/-8))
    x+=1
print("\n"+"执行结束".center(scale//2,'-'))
