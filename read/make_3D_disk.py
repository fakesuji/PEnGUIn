import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd
from matplotlib.ticker import ScalarFormatter, NullFormatter  

plt.rcParams['text.usetex'] = True

label1 = 'h27_p12J_OA_PPM4'
label2 = 'h27_p12J_a-10_OA_PPM4'

xmax = 720
ymax = 1536
zmax = 24

frame = 10

#dat1 = rd.load_3D_data('/mnt/penguin/fung/p2/',xmax,ymax,zmax,label1,frame)
dat2 = rd.load_3D_data('/mnt/penguin/fung/p2/',xmax,ymax,zmax,label2,frame)

def make_2D_density(xmax,ymax,z,den):
	den2D = np.zeros([ymax,xmax])
	for i in range(xmax):
		for j in range(ymax):
			for k in range(z.size-1):
				den2D[j,i] += (z[k+1]-z[k])*den[k,j,i]
			den2D[j,i] /= np.sqrt(np.pi/2.0)*0.027
	return den2D

#den2D1 = make_2D_density(xmax,ymax,dat1[3],dat1[4])
den2D2 = make_2D_density(xmax,ymax,dat2[3],dat2[4])
#den2D3 = make_2D_density(432,1536,dat3[3],dat3[4])
#den2D4 = make_2D_density(xmax*2,ymax*2,dat4[3],dat4[4])

#plt.figure(0)
#plt.plot(rd.cell_center(dat1[1]),np.mean(den2D1,axis=0),label="PEM3")
#plt.plot(rd.cell_center(dat2[1]),np.mean(den2D2,axis=0),label="PEM3 with OA")
#plt.plot(rd.cell_center(dat3[1]),np.mean(den2D3,axis=0),label="PEM3 1.5x")
#plt.plot(rd.cell_center(dat4[1]),np.mean(den2D4,axis=0),label="PEM3 2x")
#plt.legend()

plt.figure(0)
#plt.plot(rd.cell_center(dat1[1]),np.log10(np.mean(dat1[4][-1,:,:],axis=0)),label="PPM4")
plt.plot(rd.cell_center(dat2[1]),np.log10(np.mean(dat2[4][-1,:,:],axis=0)),label="PEM4 viscous")
plt.plot(rd.cell_center(dat2[1]),np.log10(np.mean(den2D2,axis=0)),label="PEM4 viscous surface")
plt.legend()


#plt.figure(1)
#plt.pcolormesh(dat1[1],dat1[2],np.log10(dat1[4][-1,:,:]))
#plt.colorbar()

plt.figure(2)
plt.pcolormesh(dat2[1],dat2[2],np.log10(dat2[4][-1,:,:]))
plt.colorbar()

plt.figure(3)
plt.pcolormesh(dat2[1],dat2[2],np.log10(den2D2))
plt.colorbar()

#plt.figure(3)
#plt.pcolormesh(dat3[1],dat3[2],np.log10(den2D3))
#plt.colorbar()

#plt.figure(4)
#plt.pcolormesh(dat4[1],dat4[2],np.log10(den2D4))
#plt.colorbar()

#plt.figure(5)
#tmp = dat4[4][:,ymax,:] + dat4[4][:,ymax-1,:]
#tmp/= 2
#plt.pcolormesh(dat4[1],dat4[3],np.log10(tmp))
#plt.colorbar()

#plt.figure(6)
#tmp = dat4[8][:,ymax,:] + dat4[8][:,ymax-1,:]
#tmp/= 2
#tmp/= 0.05
#zmin = np.min(tmp)
#zmax = np.max(tmp)
#norm = max(abs(zmin),abs(zmax))
#plt.pcolormesh(dat4[1],dat4[3],tmp,cmap="PRGn",vmin=-norm,vmax=norm)
#plt.colorbar()

plt.show()
