import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd
from matplotlib.ticker import ScalarFormatter, NullFormatter  

plt.rcParams['text.usetex'] = True

label = 'h25_p12J_a-10_PPM4'
xmax = 720
ymax = 1536

start = 1852

for i in range(5):
	plt.figure(i)
	dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,ymax,label,start+i)
	plt.pcolormesh(dat[1],dat[2],np.log10(dat[3]))
	plt.colorbar()

	plt.figure(5)
	plt.plot(rd.cell_center(dat[1]),np.mean(dat[3],axis=0),label='frame='+str(start+i))

plt.figure(5)
plt.legend()
plt.yscale('log')

plt.show()

#for i in range(1,1844,4):
#	dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,ymax,label,i)
#	print(i,dat[0],np.min(np.mean(dat[3],axis=0)))
	
