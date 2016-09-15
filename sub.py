import os
import random
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

def sub(a1, a2):
   os.system('./1 '+str(a1)+' '+str(a2))


os.system(' gcc SLIR_random.c -lm -o 1')

maxProcess = 30
pool = Pool(processes=maxProcess)

a1 = 0.00001
while a1 <= 0.500001:
    a2 = 0.00001
    while a2 <= 0.50001:
        pool.apply_async(sub, (a1, a2, ))
        a2 = a2 + 0.01
    a1 = a1 + 0.01

pool.close()
pool.join()
