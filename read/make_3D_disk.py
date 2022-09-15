import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd
from matplotlib.ticker import ScalarFormatter, NullFormatter  

plt.rcParams['text.usetex'] = True

label1 = 'h50_p1J_PEM3'
label2 = 'h50_p1J_OA_PEM3'

xmax = 288
ymax = 1008
zmax = 32

frame = 6

#dat1 = rd.load_3D_data('/mnt/penguin/fung/p2/',xmax,ymax,zmax,label1,frame)
dat2 = rd.load_3D_data('/mnt/penguin/fung/p2/',xmax,ymax,zmax,label2,frame)
dat3 = rd.load_3D_data('/mnt/penguin/fung/p2/',432,1536,48,label2,frame)
dat4 = rd.load_3D_data('/mnt/penguin/fung/p2/',xmax*2,ymax*2,zmax*2,label2,frame)

def make_2D_density(xmax,ymax,z,den):
	den2D = np.zeros([ymax,xmax])
	for i in range(xmax):
		for j in range(ymax):
			for k in range(z.size-1):
				den2D[j,i] += (z[k+1]-z[k])*den[k,j,i]
			den2D[j,i] /= np.sqrt(np.pi/2.0)*0.05
	return den2D

#den2D1 = make_2D_density(xmax,ymax,dat1[3],dat1[4])
den2D2 = make_2D_density(xmax,ymax,dat2[3],dat2[4])
den2D3 = make_2D_density(432,1536,dat3[3],dat3[4])
den2D4 = make_2D_density(xmax*2,ymax*2,dat4[3],dat4[4])

plt.figure(0)
#plt.plot(rd.cell_center(dat1[1]),np.mean(den2D1,axis=0),label="PEM3")
plt.plot(rd.cell_center(dat2[1]),np.mean(den2D2,axis=0),label="PEM3 with OA")
plt.plot(rd.cell_center(dat3[1]),np.mean(den2D3,axis=0),label="PEM3 1.5x")
plt.plot(rd.cell_center(dat4[1]),np.mean(den2D4,axis=0),label="PEM3 2x")
plt.legend()

#plt.figure(1)
#plt.pcolormesh(dat1[1],dat1[2],np.log10(den2D1))
#plt.colorbar()

plt.figure(2)
plt.pcolormesh(dat2[1],dat2[2],np.log10(den2D2))
plt.colorbar()

plt.figure(3)
plt.pcolormesh(dat3[1],dat3[2],np.log10(den2D3))
plt.colorbar()

plt.figure(4)
plt.pcolormesh(dat4[1],dat4[2],np.log10(den2D4))
plt.colorbar()

plt.figure(5)
tmp = dat4[4][:,ymax,:] + dat4[4][:,ymax-1,:]
tmp/= 2
plt.pcolormesh(dat4[1],dat4[3],np.log10(tmp))
plt.colorbar()

plt.figure(6)
tmp = dat4[8][:,ymax,:] + dat4[8][:,ymax-1,:]
tmp/= 2
tmp/= 0.05
zmin = np.min(tmp)
zmax = np.max(tmp)
norm = max(abs(zmin),abs(zmax))
plt.pcolormesh(dat4[1],dat4[3],tmp,cmap="PRGn",vmin=-norm,vmax=norm)
plt.colorbar()

plt.show()
