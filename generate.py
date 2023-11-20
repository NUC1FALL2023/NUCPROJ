import numpy as np
from tabulate import tabulate
import time
from HartreeFock import *
from CGgenerator import *

    
##################################################################
###                  Generate Basis States                     ###
##################################################################

# Array of parameters for gradient descent
x = []

ite = 0 # counting dummy variable

# creating a file for any given 'n' value
f = open("sdbasisGeneratedRaw.dat", 'w')
data = []
eps = []
statesDict = []

# States that should be generated
# n = [0,1,0]
# l = [2,0,2]
# j = [5,1,3]

# Start with lowest lying states
n = [0,1]
l = [2,0]
j = [5,1]

nlist = []
llist = []
jlist = []
mjlist = []
            
# Optimized Energies from literature 
# realEps = -np.asarray([4.15,4.15,8.02,8.02,3.94,3.94,7.63,7.63,3.77,3.77,6.84,6.84,2.75,2.75,3.61,3.61,8.62,8.62,11.42,11.42,7,7,9.9,9.9])
realEps = -np.asarray([4.15,4.15,4.15,4.15,4.15,4.15,2.75,2.75])

# next we need to go over all possible combinations of quantum numbers, starting with 'n'
for i in range(len(n)): # there is only l=0 and l=2 in this case but could be adjusted for any number
        mj_min = -j[i] # mj goes from j to -j in increments of 1
        mj = j[i]
        while mj >= mj_min:
            # val = np.random.uniform(-2,2)
            val = -realEps[ite]
            data.append([ite,n[i],l[i],j[i],int(mj),1,val])
            eps.append(val)
            ite +=1
            # data.append([ite,n[i],l[i],j[i],int(mj),-1,val])
            # ite +=1
            statesDict.append(dict(n=int(n[i]),l=int(l[i]),j=int(j[i]),mj=int(mj),t3=int(1)))
            # statesDict.append(dict(n=int(n[i]),l=int(l[i]),j=int(j[i]),mj=int(mj),t3=int(-1)))
            nlist.append(n[i])
            llist.append(l[i])
            jlist.append(j[i])
            mjlist.append(mj)
            mj = mj-2
            # eps.append(val)
            # x.append(val)

# the tabulate function is slow but that should not matter for these calculations
f.write(tabulate(data, tablefmt="plain",showindex=False))
f.close()

converge = [7.767,7.767,7.767,7.767,7.767,7.767,7.767,7.767]
    
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
                    
                        if a < j[0]+1 and b < j[0]+1 and c < j[0]+1 and d < j[0]+1:
                            ME = np.random.normal(3.91,3.91*0.12,1)
                        elif a > j[0] and b > j[0] and c > j[0] and d > j[0]:
                            ME = np.random.normal(3.65,3.65*0.12,1)
                        else:
                            ME = np.random.normal(2.68,2.68*0.12,1)
                            
                        CG = 0
                        Mp = int(mjlist[c]/2+mjlist[d]/2)
                        Jp = int(jlist[c]/2+jlist[d]/2)
                        for JM in range(Mp,Jp+1):
                            CG += CFcoeff(jlist[c],jlist[d],mjlist[c],mjlist[d],JM,Mp)
                            
                        ME = CG*ME
                        x.append(ME)
                        
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
 
#print(end-start)

    
##################################################################
###        Optimize Energies and Two Body Matrix Elements      ###
##################################################################

step = 1e-3
# rate = 1e-3
rate = np.zeros(len(x))
A = 2
Error = 1e20
dE = 1e-5
stepN = 0
maxSteps = 1000
max_iter = 20

def getRate(N):
    return np.floor(np.log10(abs(N)))

def x2epsV(vec_x):
    energy = np.zeros(ite)
    tbme = np.zeros([ite,ite,ite,ite])
    
    ind = 0
    
# =============================================================================
#     for i in range(int(ite/2)):
#         energy[2*i] = vec_x[ind]
#         energy[2*i+1] = vec_x[ind]
#         
#         ind += 1
# =============================================================================
        
    for a in range(ite):
        for b in range(ite):
            if (a < b):
                for c in range(ite):
                    for d in range(ite):
                        if (c < d):
                            
                            ME = vec_x[ind]
                            ind += 1
                                
                            tbme[a][b][c][d] = ME
                            tbme[b][a][d][c] = ME
                            tbme[b][a][c][d] = -ME
                            tbme[a][b][d][c] = -ME
                            
    return energy, tbme




def calcError(E):
    return np.sum(np.square(E[:A]-converge[:A]))


oldError = 0
diffError = 1
        
# Find the gradient for each element

while diffError > dE and stepN < maxSteps:

    df = np.zeros(len(x))

    for i in range(len(x)):
        
        xPlus = x.copy()
        xMinus = x.copy()
        
        xPlus[i] += step
        xMinus[i] -= step
        
        skip, VPlus = x2epsV(xPlus)
        skip, VMinus = x2epsV(xMinus)
        
        epsNewPlus, skip = HartreeFock(A,1e-5,max_iter,statesDict,eps,VPlus)    
        epsNewMinus, skip = HartreeFock(A,1e-5,max_iter,statesDict,eps,VMinus)
        
        fplus = calcError(epsNewPlus[-1])
        fminus = calcError(epsNewMinus[-1])
        
        df[i] = (fplus-fminus)/(2*step)
        
        rate[i] = 10**(-(getRate(df[i])+3))
        
    x = x-rate*df
    
    skip, VN = x2epsV(x)
    
    epsNew, skip = HartreeFock(A,1e-5,max_iter,statesDict,eps,VN)
    
    oldError = Error
    Error = calcError(epsNew[-1])
    
    diffError = abs(Error-oldError)
    
    print(Error)
    print(epsNew[-1][:A])
    
    stepN += 1
    
skip, finalV = x2epsV(x)

finalEps = epsNew.copy()

f = open("sdbasisGeneratedOptimized.dat","w")
dataOpt = []
k = 0

for i in range(len(n)): # there is only l=0 and l=2 in this case but could be adjusted for any number
        mj_min = -j[i] # mj goes from j to -j in increments of 1
        mj = j[i]
        while mj >= mj_min:
            dataOpt.append([ite,n[i],l[i],j[i],int(mj),1,finalEps[k]])
            k += 1
            dataOpt.append([ite,n[i],l[i],j[i],int(mj),-1,finalEps[k]])
            k += 1
            mj = mj-2
            

f.write(tabulate(dataOpt, tablefmt="plain",showindex=False))
f.close()
    
f = open("tbmeGeneratedOptimized.dat","w")

for a in range(ite):
    for b in range(ite):
        if (a < b):
            for c in range(ite):
                for d in range(ite):
                    if (c < d):
                        
                        f.write(str(a) + '\t' + str(b) + '\t'+ str(c) + '\t' + str(d) + '\t' + str(finalV[a][b][c][d]) + '\n')
                        f.write(str(b) + '\t' + str(a) + '\t'+ str(c) + '\t' + str(d) + '\t' + str(finalV[b][a][c][d]) + '\n')
                        f.write(str(a) + '\t' + str(b) + '\t'+ str(d) + '\t' + str(c) + '\t' + str(finalV[a][b][d][c]) + '\n')
                        f.write(str(b) + '\t' + str(a) + '\t'+ str(d) + '\t' + str(c) + '\t' + str(finalV[b][a][d][c]) + '\n')
                        
f.close()