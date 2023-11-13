import numpy as np
from tabulate import tabulate
import time

    
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
j = [3,1,5]


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
            
# Optimized Energies from literature (placeholder until we find these)
realEps = np.zeros(ite)

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

    
##################################################################
###                   Optimize Energies                        ###
##################################################################

dE = 1e-5
step = 0.0001
max_iter = 5

minError = 1e20
newEps = eps.copy()

for i in range(ite):
    
    while minError > dE:

        epsPlus = eps.copy()
        epsMinus = eps.copy()
        epsPlus[i] += step
        epsMinus[i] -= step
        
        bindingEnergiesStart, skip = HartreeFock(16,1e-5,max_iter,statesDict,eps,V)    
        bindingEnergiesPlus, skip = HartreeFock(16,1e-5,max_iter,statesDict,epsPlus,V)    
        bindingEnergiesMinus, skip = HartreeFock(16,1e-5,max_iter,statesDict,epsMinus,V)
        
        dEStart = np.average(abs(bindingEnergiesStart - realEps))
        dEPlus = np.average(abs(bindingEnergiesPlus - realEps))
        dEMinus = np.average(abs(bindingEnergiesMinus - realEps))
    
        if dEPlus < dEStart:
            if dEMinus < dEStart:
                newEps = epsMinus.copy()
                minError = dEMinus
            else:
                newEps = epsPlus.copy()
                minError = dEPlus
        else: 
            break
        
        eps = newEps.copy() 
    
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
dE = 1e-5
step = 0.0001
max_iter = 5

minError = 1e20
newV = V.copy()

VPlus = V.copy()
VMinus = V.copy()

for a in range(ite):
    for b in range(ite):
        if (a < b):
            for c in range(ite):
                for d in range(ite):
                    if (c < d):
                        while minError > dE:

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
        
                            bindingEnergiesStart, skip = HartreeFock(16,1e-5,max_iter,statesDict,newEps,V)    
                            bindingEnergiesPlus, skip = HartreeFock(16,1e-5,max_iter,statesDict,newEps,VPlus)    
                            bindingEnergiesMinus, skip = HartreeFock(16,1e-5,max_iter,statesDict,newEps,VMinus)
        
                            dEStart = np.average(abs(bindingEnergiesStart - realEps))
                            dEPlus = np.average(abs(bindingEnergiesPlus - realEps))
                            dEMinus = np.average(abs(bindingEnergiesMinus - realEps))
                            
                            if dEPlus < dEStart:
                                if dEMinus < dEStart:
                                    newV = VMinus.copy()
                                    minError = dEMinus
                                else:
                                    newV = VPlus.copy()
                                    minError = dEPlus
                            else: 
                                break
                            
                            V = newV.copy()
                            
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