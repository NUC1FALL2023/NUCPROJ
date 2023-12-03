import numpy as np
from HartreeFock import *

##################################################################
###             READ DAT FILES PROVIDED                        ###
##################################################################

# Data is stored in dictionaries for explicit labeling

n = [0,1]
l = [2,0]
j = [5,1]

statesDict = []

realEps = np.asarray([4.15,4.15,4.15,4.15,4.15,4.15,2.75,2.75])
for i in range(len(n)): # there is only l=0 and l=2 in this case but could be adjusted for any number
        mj_min = -j[i] # mj goes from j to -j in increments of 1
        mj = j[i]
        while mj >= mj_min:
            statesDict.append(dict(n=int(n[i]),l=int(l[i]),j=int(j[i]),mj=int(mj),t3=int(1)))
            mj = mj-2
    
lenStates = len(statesDict)
    
Vfile = open('tbmeGeneratedRaw.dat','r')

VDict = []

for x in Vfile:
    temp1,temp2,temp3,temp4,temp5 = x.split()
    if temp5 != "nan":
        VDict.append(dict(a=int(temp1),b=int(temp2),c=int(temp3),d=int(temp4),val=float(temp5)))
    
##################################################################
###                 INITIALIZING MATRICES                      ###
##################################################################
    
# Setting up matrices
    
V = np.zeros((lenStates,lenStates,lenStates,lenStates))

for x in VDict:
    V[x["a"]][x["b"]][x["c"]][x["d"]] = x["val"]

eps = realEps
    
##################################################################
###                    Run Hartree Fock                        ###
##################################################################
    
# Running Hartree Fock
# Function is stored in HartreeFock.py
# Inputs are HartreeFock( A , Smallest Difference in Energy Wanted ,
#                           Max Number of iterations requests , 
#                           Basis Stored as a Dict , Energies in an Array,
#                           and Two Body Density Matrix as 4d Array)
# Outputs are bindingEnergies, a 2d array where the first index is the iteration number
#                           and the second index is the state index
# To access the final iteration, use bindingEnergies[-1] as the array
# Expected negEnergies = A

# A is number of nucleons
A = np.asarray([1,2,3,4,5,6])

for AN in A:
    bindingEnergies, negEnergies = HartreeFock(AN,1e-5,20,statesDict,eps,V)
    
    # Writing out final energies
    # orbit defined by l quantum number
    orbit = ["s","p","d","f","g","h"]
    f = open("BindingEnergiesHO.txt","w")
    
    for i in range(lenStates):
        # prints out in orbital notation nl^pi,nu (j)
        if statesDict[i]["t3"] == 1: nuc = "n"
        else: nuc = "p"
        f.write(str(statesDict[i]["n"])+orbit[statesDict[i]["l"]]+"^"+nuc+"("+str(statesDict[i]["j"])+"/2) " + str(round(bindingEnergies[-1][i],4))+ " MeV\n")
        
    B = np.zeros(len(bindingEnergies))
    
    # Calculates the binding energy per nucleon and stores it in an array
    for i in range(len(bindingEnergies)):
        B[i] = sum(bindingEnergies[i][:negEnergies])/negEnergies
            
    print((np.average(bindingEnergies[-1][:AN])*AN-7.976*16)/(16+AN))
    
    f.close()
