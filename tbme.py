import numpy as np
from tabulate import tabulate
import time

start = time.time()
q = int(input("Enter quanta number: "))
f = open("tbme_" + str(q)+ ".txt", 'w')
table = []

for a in range(1,q+1):
	for b in range(1,q+1):
		if (a != b and a < b):
			for c in range(1,q+1):
				for d in range(1, q+1):
					if (c !=d and c < d):
						val = np.random.uniform(-2,2)
						f.write(str(a) + '\t' + str(b) + '\t'+ str(c) + '\t' + str(d) + '\t' + str(val) + '\n')
						f.write(str(b) + '\t' + str(a) + '\t'+ str(c) + '\t' + str(d) + '\t' + str(-val) + '\n')
						f.write(str(a) + '\t' + str(b) + '\t'+ str(d) + '\t' + str(c) + '\t' + str(-val) + '\n')
						f.write(str(b) + '\t' + str(a) + '\t'+ str(d) + '\t' + str(c) + '\t' + str(val) + '\n')
						#table.append([a,b,c,d,val])
						#table.append([b,a,c,d,-val])
						#table.append([a,b,d,c,-val])
						#table.append([b,a,d,c,val])
#print(len(table))
#f.write(tabulate(table, tablefmt="plain"))
f.close()
end = time.time()

print(end-start)
