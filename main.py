##################################################################
###             READ DAT FILES PROVIDED                        ###
##################################################################
# Data is stored in dictionaries for explicit labeling

statefile = open('spbasis.dat','r')

states = []

for x in statefile:
    temp1,temp2,temp3,temp4,temp5,temp6 = x.split()
    states.append(dict(n=temp2,l=temp3,j=temp4,mj=temp5,t3=temp6))
    
Vfile = open('tbme.dat','r')

V = []

for x in Vfile:
    temp1,temp2,temp3,temp4,temp5 = x.split()
    V.append(dict(a=temp1,b=temp2,c=temp3,d=temp4,val=temp5))


##################################################################
###                         END                                ###
##################################################################
