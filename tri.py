import numpy as np
import matplotlib.pyplot as plt
import copy
import math
from pprint import pprint
from scipy.ndimage import convolve, generate_binary_structure

energy1 = []; energy2=[];energy3=[];energy4=[];energy5=[];

energy = [energy1,energy2,energy3,energy4,energy5]

energytxt = ["E_MN_1000.txt","E_MN_10000.txt","E_MN_100000.txt","E_MN_1000000.txt","E_MN_10000000.txt"]


for i in range(len(energy)):
	with open(energytxt[i],"r") as f:
		for line in f:
			energy[i].append(float(line.strip()))

plt.figure(figsize=(12.8,9.6))
plt.title('Energy(E)vs Temperature(T)')

x = np.arange(0.1, 5.1, 0.1)

plt.xticks(np.arange(0.1,5.5,0.5))
plt.xlabel('Temperature(T)')
plt.ylabel('Energy(E)')

for i in range(len(energy)):
	y = energy[i]
	plt.scatter(x,y)
	plt.plot(x,y,label='MN=%d'%(1000*10**i))

plt.legend(loc='lower right')
plt.show()
