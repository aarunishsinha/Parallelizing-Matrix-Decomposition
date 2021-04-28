import numpy as np
import sys

n = int(sys.argv[1])
file = f'input_{n}_{n}.txt'

A = np.random.uniform(-4, 4, (n, n))

f = open(file, 'w')
for i in A:
    for j in i:
        print ("{0:.12f}".format(j), '', end='', file=f)
    print(file=f)
f.close()