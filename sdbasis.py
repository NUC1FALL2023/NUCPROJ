import numpy as np
from tabulate import tabulate
import pandas as pd

n_max = int(input('Give n(max like 3 tbh): '))
l = [0,2]
#s = 1
#j = l + s
#mj = -j to j
#t = +/- 1 for each state
ite = 1
f = open("sdbasis_nmax" + str(n_max) + ".txt", 'w')
data = []

for i in range(len(l)):
	n = 0
	while n <= n_max:
		j = l[i] + 1 / 2
		mj_min = -j
		mj = j
		while mj >= mj_min:
			for c in range(2):
				if c == 0:
					data.append([ite,n,l[i],int(l[i]+ 1 % 2),int(2*mj),1])
					ite +=1
				if c == 1:
					data.append([ite,n,l[i],int(l[i] + 1 % 2),int(2* mj),-1])
					ite +=1
			mj = mj-1
		n += 1	

f.write(tabulate(data, tablefmt="plain",showindex=False))
f.close()

