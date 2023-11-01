##################################################################
###             READ DAT FILES PROVIDED                        ###
##################################################################
# Data is stored in dictionaries for explicit labeling

statefile = open('spbasis.dat','r')

states = []

for x in statefile:
    temp1,temp2,temp3,temp4,temp5,temp6 = x.split()
    states.append(dict(n=int(temp2),l=int(temp3),j=int(temp4),mj=int(temp5),t3=int(temp6)))
    
Vfile = open('tbme.dat','r')

V = []

for x in Vfile:
    temp1,temp2,temp3,temp4,temp5 = x.split()
    V.append(dict(a=int(temp1),b=int(temp2),c=int(temp3),d=int(temp4),val=float(temp5)))
    
for s in states:
    s.update({"N":2*s["n"]+s["l"]})

##################################################################
###                         END                                ###
##################################################################
