#! /bin/python
import os

file = open("temp.list")

filelist=[]

while 1:
    lines = file.readlines(100000)
    if not lines:
        break
    for line in lines:
        linne=line[-32:-1]
        print(linne)
        filelist.append(linne)

        pass # do something



tuplelink=os.walk("/Users/yuanjing/Work/STAR/Hypertron/H3L_3body_yuehang/data/mixed/out_KF")
for tuple in tuplelink:
    pass
listlink=tuple[2]
misslist=[]

l=len(listlink)
for i in range(l):
    try:
        print(listlink.index(filelist[i]),filelist[i],'exists')
    except:
        print(listlink[i],'has missed')
        misslist.append(listlink[i])

print(misslist,'are missed')

