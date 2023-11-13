import numpy as np
from tabulate import tabulate
import time
from HartreeFock import *

    
##################################################################
###                  Generate Basis States                     ###
##################################################################

ite = 0 # counting dummy variable

# creating a file for any given 'n' value
f = open("sdbasisGeneratedRaw.dat", 'w')
data = []
eps = []
statesDict = []

# States that should be generated
n = [0,1,0]
l = [2,0,2]
j = [5,1,3]


# next we need to go over all possible combinations of quantum numbers, starting with 'n'
for i in range(len(n)): # there is only l=0 and l=2 in this case but could be adjusted for any number
        mj_min = -j[i] # mj goes from j to -j in increments of 1
        mj = j[i]
        while mj >= mj_min:
            val = np.random.uniform(-2,2)
            data.append([ite,n[i],l[i],j[i],int(mj),1,val])
            eps.append(val)
            ite +=1
            data.append([ite,n[i],l[i],j[i],int(mj),-1,val])
            ite +=1
            statesDict.append(dict(n=int(n[i]),l=int(l[i]),j=int(j[i]),mj=int(mj),t3=int(1)))
            statesDict.append(dict(n=int(n[i]),l=int(l[i]),j=int(j[i]),mj=int(mj),t3=int(-1)))
            mj = mj-2
            eps.append(val)
            
# Optimized Energies from literature 
realEps = [4.15,4.15,8.02,8.02,3.94,3.94,7.63,7.63,3.77,3.77,6.84,6.84,2.75,2.75,3.61,3.61,8.62,8.62,11.42,11.42,7,7,9.9,9.9]

# the tabulate function is slow but that should not matter for these calculations
f.write(tabulate(data, tablefmt="plain",showindex=False))
f.close()

    
##################################################################
###           Generate Two Body Matrix Elements                ###
##################################################################


start = time.time() # I added this to help optimize the code
 
f = open("tbmeGeneratedRaw.dat", 'w')

V = np.random.rand(ite,ite,ite,ite)*2-1

for a in range(ite):
    for c in range(ite):
        for d in range(ite):
            # Antisymmetry forces these "diagonal" terms to be 0
            V[a][a][c][d] = 0
            V[a][a][d][c] = 0
            V[c][d][a][a] = 0
            V[d][c][a][a] = 0
            
for a in range(ite):
    for b in range(ite):
        if (a < b):
            for c in range(ite):
                for d in range(ite):
                    if (c < d):
                        # Normalized and antisymmetric
                        ME = 1/2*(V[a][b][c][d]-V[a][b][d][c]-V[b][a][c][d]+V[b][a][d][c])
                        V[a][b][c][d] = ME
                        V[b][a][d][c] = ME
                        V[b][a][c][d] = -ME
                        V[a][b][d][c] = -ME
                        
                        f.write(str(a) + '\t' + str(b) + '\t'+ str(c) + '\t' + str(d) + '\t' + str(ME) + '\n')
                        f.write(str(b) + '\t' + str(a) + '\t'+ str(c) + '\t' + str(d) + '\t' + str(-ME) + '\n')
                        f.write(str(a) + '\t' + str(b) + '\t'+ str(d) + '\t' + str(c) + '\t' + str(-ME) + '\n')
                        f.write(str(b) + '\t' + str(a) + '\t'+ str(d) + '\t' + str(c) + '\t' + str(ME) + '\n')
 						
f.close()
end = time.time()
 
print(end-start)

# Options for optimization
dE = 1e-5
step = 0.0001
max_iter = 5
A = 1
max_steps = 100
    
##################################################################
###                   Optimize Energies                        ###
##################################################################

Error = 1e20
newEps = eps.copy()

for i in range(ite):
    
    step_iter = 0
    
    while Error > dE and step_iter < max_steps:

        epsPlus = eps.copy()
        epsMinus = eps.copy()
        epsPlus[i] += step
        epsMinus[i] -= step
        
        bindingEnergiesStart, skip = HartreeFock(A,1e-5,max_iter,statesDict,eps,V)    
        bindingEnergiesPlus, skip = HartreeFock(A,1e-5,max_iter,statesDict,epsPlus,V)    
        bindingEnergiesMinus, skip = HartreeFock(A,1e-5,max_iter,statesDict,epsMinus,V)
        
        dEStart = np.sum(np.square(bindingEnergiesStart[-1] - realEps))
        dEPlus = np.sum(np.square(bindingEnergiesPlus[-1] - realEps))
        dEMinus = np.sum(np.square(bindingEnergiesMinus[-1] - realEps))
    
        fp = (dEPlus-dEMinus)/(2*step)
        f = dEStart
        
        newEps = eps.copy()
        newEps[i] = eps[i] - fp/f
        
        eps = newEps.copy()
        
        checkBindingEnergy, skip = HartreeFock(A,1e-5,max_iter,statesDict,eps,V)
        Error = np.sum(np.square(checkBindingEnergy[-1] - realEps))
        step_iter += 1
        
    print("Converged for energy " + str(i) + "\n")
    
f = open("sdbasisGeneratedOptimized.dat","w")
dataOpt = []
k = 0

for i in range(len(n)): # there is only l=0 and l=2 in this case but could be adjusted for any number
        mj_min = -j[i] # mj goes from j to -j in increments of 1
        mj = j[i]
        while mj >= mj_min:
            dataOpt.append([ite,n[i],l[i],j[i],int(mj),1,newEps[k]])
            k += 1
            dataOpt.append([ite,n[i],l[i],j[i],int(mj),-1,newEps[k]])
            k += 1
            mj = mj-2
            

f.write(tabulate(dataOpt, tablefmt="plain",showindex=False))
f.close()


    
##################################################################
###          Optimize Two Body Matrix Elements                 ###
##################################################################

Error = 1e20
newV = V.copy()

VPlus = V.copy()
VMinus = V.copy()

for a in range(ite):
    for b in range(ite):
        if (a < b):
            for c in range(ite):
                for d in range(ite):
                    if (c < d):
                        
                        step_iter = 0
                        
                        while Error > dE and step_iter < max_steps:

                            VPlus = V.copy()
                            VMinus = V.copy()
                            
                            VPlus[a][b][c][d] += step
                            VPlus[b][a][d][c] += step
                            VPlus[b][a][c][d] -= step
                            VPlus[a][b][d][c] -= step
                            
                            VMinus[a][b][c][d] -= step
                            VMinus[b][a][d][c] -= step
                            VMinus[b][a][c][d] += step
                            VMinus[a][b][d][c] += step
        
                            bindingEnergiesStart, skip = HartreeFock(A,1e-5,max_iter,statesDict,newEps,V)    
                            bindingEnergiesPlus, skip = HartreeFock(A,1e-5,max_iter,statesDict,newEps,VPlus)    
                            bindingEnergiesMinus, skip = HartreeFock(A,1e-5,max_iter,statesDict,newEps,VMinus)
        
                            dEStart = np.sum(np.square(bindingEnergiesStart[-1] - realEps))
                            dEPlus = np.sum(np.square(bindingEnergiesPlus[-1] - realEps))
                            dEMinus = np.sum(np.square(bindingEnergiesMinus[-1] - realEps))
                            
                            fp = (dEPlus-dEMinus)/(2*step)
                            f = dEStart
                            
                            newV = V.copy()
                            newEps[a][b][c][d] = V[a][b][c][d] - fp/f
                            newEps[b][a][d][c] = V[b][a][d][c] - fp/f
                            newEps[a][b][d][c] = V[a][b][d][c] + fp/f
                            newEps[b][s][c][d] = V[b][s][c][d] + fp/f
                            
                            V = newV.copy()
                            
                            checkBindingEnergy, skip = HartreeFock(A,1e-5,max_iter,statesDict,newEps,V)
                            Error = np.sum(np.square(checkBindingEnergy[-1] - realEps))
                            step_iter += 1
                            
f = open("tbmeGeneratedOptimized.dat","w")

for a in range(ite):
    for b in range(ite):
        if (a < b):
            for c in range(ite):
                for d in range(ite):
                    if (c < d):
                        
                        f.write(str(a) + '\t' + str(b) + '\t'+ str(c) + '\t' + str(d) + '\t' + str(newV[a][b][c][d]) + '\n')
                        f.write(str(b) + '\t' + str(a) + '\t'+ str(c) + '\t' + str(d) + '\t' + str(newV[b][a][c][d]) + '\n')
                        f.write(str(a) + '\t' + str(b) + '\t'+ str(d) + '\t' + str(c) + '\t' + str(newV[a][b][d][c]) + '\n')
                        f.write(str(b) + '\t' + str(a) + '\t'+ str(d) + '\t' + str(c) + '\t' + str(newV[b][a][d][c]) + '\n')
                        
f.close()