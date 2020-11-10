import matplotlib.pyplot as plt
import numpy as np

x=[1,2,3,4,5,6]
y=[3,9,1,6,4,7]
z=[8,3,6,2,5,7]
#t=[1,2,3,4,5,6,7,8,9,10]

plt.plot(x,y) #常规折线图
plt.xlabel('xlabel')
plt.ylabel('ylabel')
plt.show()

plt.plot(x,z,'ro') #散点图,'ro'表示红色的圆点
plt.axis([0, 10, 0, 10]) #设置x与y轴的长度
plt.show()

t = np.arange(0.0, 5.0, 0.2) #从0到5步长0.2生成序列
plt.plot(t, t, 'r--', t, t**2, 'bs', t, t**3, 'g^') # 复合散点图,'r--':红色的需要;'bs':蓝色方块;'g^':绿色三角
plt.show()
