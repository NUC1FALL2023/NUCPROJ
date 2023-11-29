import numpy as np
from math import factorial as fac

def CGcoeff(j1,j2,m1in,m2in,J,Min):
    
    if Min<0:
        M = -Min
        m1 = -m1in
        m2 = -m2in
    else:
        M = Min
        m1 = m1in
        m2 = m2in
    
    CG = 0
    
    if M == m1+m2:
        if J+j1-j2 < 0 or J-j1+j2 < 0 or j1+j2-J < 0:
           return 0
        CG += np.sqrt((int(2*J+1))*fac(int(J+j1-j2))*fac(int(J-j1+j2))*fac(int(j1+j2-J))/fac(int(j1+j2+J+1)))
        
    CG *= np.sqrt(fac(int(J+M))*fac(int(J-M))*fac(int(j1-m1))*fac(int(j1+m1))*fac(int(j2-m2))*fac(int(j2+m2)))
    
    k1 = j1+j2-J
    k2 = j1-m1
    k3 = j2+m2
    
    if k1 <= k2 and k1 <= k3:
        kmax = k1
    elif k2 <= k1 and k2 <= k3:
        kmax = k2
    else:
        kmax = k3     
    
    sumk = 0
    
    for k in range(int(kmax)+1):
        if J-j2+m1+k >= 0 and J-j1-m2+k >= 0:
            sumk += np.power(-1,k)/(fac(k)*fac(int(j1+j2-J-k))*fac(int(j2+m2-k))*fac(int(J-j2+m1+k))*fac(int(J-j1-m2+k)))
    
    finalCG = CG*sumk
    
    if Min<0:
        finalCG = finalCG*np.power(-1,J-j1-j2)
    
    return finalCG