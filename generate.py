import numpy as np
from tabulate import tabulate
import time

ite = 0 # counting dummy variable

# creating a file for any given 'n' value
f = open("sdbasisGenerated.dat", 'w')
data = []

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
            ite +=1
            data.append([ite,n[i],l[i],j[i],int(mj),-1,val])
            ite +=1
            mj = mj-2

# the tabulate function is slow but that should not matter for these calculations
f.write(tabulate(data, tablefmt="plain",showindex=False))
f.close()

start = time.time() # I added this to help optimize the code
 
f = open("tbmeGenerated.dat", 'w')

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

# Optimize Energies and Two-Body Matrix Elements

