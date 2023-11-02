import numpy as np

##################################################################
###             READ DAT FILES PROVIDED                        ###
##################################################################

# Data is stored in dictionaries for explicit labeling

statefile = open('spbasis.dat','r')

statesDict = []

for x in statefile:
    temp1,temp2,temp3,temp4,temp5,temp6 = x.split()
    statesDict.append(dict(n=int(temp2),l=int(temp3),j=int(temp4),mj=int(temp5),t3=int(temp6)))
    
lenStates = len(statesDict)
    
Vfile = open('tbme.dat','r')

VDict = []

for x in Vfile:
    temp1,temp2,temp3,temp4,temp5 = x.split()
    VDict.append(dict(a=int(temp1)-1,b=int(temp2)-1,c=int(temp3)-1,d=int(temp4)-1,val=float(temp5)))
    
for s in statesDict:
    s.update({"N":2*s["n"]+s["l"]})    
    
##################################################################
###                 INITIALIZING MATRICES                      ###
##################################################################
    
# Setting up matrices
    
V = np.zeros((lenStates,lenStates,lenStates,lenStates))

for x in VDict:
    V[x["a"]][x["b"]][x["c"]][x["a"]] = x["val"]

# Energies of Harmonic Oscillator
eps = np.zeros(lenStates)

hw = 1 # HO quanta, chosen in MeV
for i in range(lenStates):
    eps[i] = hw*(statesDict[i]["N"]+3/2)
    
# Take a random matrix to guess C[delta][gamma]
C = np.random.rand(lenStates,lenStates)

##################################################################
##################################################################
##################################################################
###                 START OF THE HF METHOD                     ###
##################################################################
##################################################################
##################################################################

# Number of nucleons
A = 16

# How small of a difference in energy do we want?
saturation = 0.001 
dE = 1
max_iter = 10
it_count = 0
ener_old = np.zeros(lenStates)
energies = np.zeros(lenStates)

while ((dE > saturation) and (it_count < max_iter)):
    # Calculate the one body density matrix
    rho = np.zeros((lenStates,lenStates))

    for delta in range(lenStates):
        for gamma in range(lenStates):
            for j in range(A):
                rho[delta][gamma] += C[delta][j]*C[gamma][j]
    
    # Calculate HF Potential
    U = np.zeros((lenStates,lenStates))
    
    for eta in range(lenStates):
        for beta in range(lenStates):
            for delta in range(lenStates):
                for gamma in range(lenStates):
                    U[eta][beta] += rho[delta][gamma]*V[eta][gamma][beta][delta]
                    
    # For initial <alpha|h_0|beta>, use diagonal matrix with energies
    E = np.diag(eps)
            
    # Calculate Hartree-Fock Hamiltonian
    HF = E+U
    
    # Calculate new C matrix and energies, store previous energies
    ener_old = energies
    energies, C = np.linalg.eigh(HF)
    
    # Calculate energy difference
    dE = np.sum(np.abs(energies-ener_old))
    it_count += 1
    
    print("Iteration " + str(it_count))
    for i in range(lenStates):
        print(str(statesDict[i]["n"]) + " " + str(statesDict[i]["l"]) + " " + str(statesDict[i]["j"]) + " " + str(statesDict[i]["mj"]) + " " + str(statesDict[i]["t3"]) + " " + str(energies[i]))
        
    print("\n")
