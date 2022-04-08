import read as pg
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

d1 = pg.load_2D_data("/mnt/penguin/fung/p2/",960,2016,"h25_p12J_a-1_OA_PLM",0)

for i in range(0,1):
	m1 = pg.load_2D_data("/mnt/penguin/fung/p2/",960,2016,"h25_p12J_a-1_OA_PLM",i)
	print(np.max(np.abs(d1[3]-m1[3])))
	plt.figure(i)
	#plt.pcolormesh(d1[1],d1[2],m1[3],norm=matplotlib.colors.LogNorm())
	plt.pcolormesh(d1[1],d1[2],m1[3],vmax=10.0,vmin=0.0)
	plt.colorbar()
	plt.savefig("../images/CItau_"+str(i)+".png")
	plt.close()
