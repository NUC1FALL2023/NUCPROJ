##################################################################
###             READ DAT FILES PROVIDED                        ###
##################################################################
# Data is stored in dictionaries for explicit labeling

statefile = open('spbasis.dat','r')

n = []
l = []
j = []
m = []
t = []

for x in statefile:
    temp1,temp2,temp3,temp4,temp5,temp6 = x.split()
    n.append(int(temp2))
    l.append(int(temp3))
    j.append(int(temp4))
    m.append(int(temp5))
    t.append(int(temp6))
    
Vfile = open('tbme.dat','r')

a = []
b = []
c = []
d = []
Vabcd = []

for x in Vfile:
    temp1,temp2,temp3,temp4,temp5 = x.split()
    a.append(int(temp1))
    b.append(int(temp2))
    c.append(int(temp3))
    d.append(int(temp4))
    Vabcd.append(float(temp5))
    
states = []
for i in range(len(n)):
    states.append(dict(n=n[i],l=l[i],j=j[i],mj=m[i],t3=t[i]))
    
V = []
for i in range(len(Vabcd)):
    V.append(dict(a=a[i],b=b[i],c=c[i],d=d[i],val=Vabcd[i]))

##################################################################
###                         END                                ###
##################################################################
