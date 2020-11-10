temp=input()
if temp[0] in ['C']:
    ans=1.8*eval(temp[1:-1]+temp[-1])+32  #TempStr[1:] 除第一位（0号位）之外的所有字符
    print("F{:.2f}".format(ans))
elif temp[0] in ['F']:
    ans=(eval(temp[1:-1]+temp[-1])-32)/1.8
    print("C{:.2f}".format(ans))
else:
    print('输入格式错误')
