import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd


label = 'h35_p2J_PPM'
xmax = 192
ymax = 672
zmax = 72

frame = 7
dat0 = rd.load_3D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,zmax,label,0)
p0 = rd.load_3D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,zmax,label,frame)

dr = dat0[1][1:xmax+1]-dat0[1][0:xmax]
dp = dat0[2][1:ymax+1]-dat0[2][0:ymax]
dz = dat0[3][1:zmax+1]-dat0[3][0:zmax]

diff = dat0[4]

for i in range(0,xmax):
	for j in range(0,ymax):
		diff[:,j,i] = p0[4][:,j,i]*dz
plt.figure(1)
plt.title("vertically integrated density")
plt.pcolormesh(dat0[1],dat0[2],np.log10(np.sum(diff,axis=0)/0.035/np.sqrt(np.pi/2.0)))
plt.colorbar()

for i in range(0,xmax):
	for j in range(0,ymax):
		diff[:,j,i] = p0[5][:,j,i]*dz
plt.figure(2)
plt.title("vertically integrated pressure")
plt.pcolormesh(dat0[1],dat0[2],np.log10(np.sum(diff,axis=0)/0.035/np.sqrt(np.pi/2.0)))
plt.colorbar()

for i in range(0,xmax):
	for j in range(0,ymax):
		diff[:,j,i] = p0[6][:,j,i]*dz
plt.figure(3)
plt.title("vertically averaged vr")
plt.pcolormesh(dat0[1],dat0[2],np.sum(diff,axis=0)/np.sum(dz))
plt.colorbar()

for i in range(0,xmax):
	for j in range(0,ymax):
		diff[:,j,i] = p0[7][:,j,i]*dz
plt.figure(4)
plt.title("vertically averaged vphi")
plt.pcolormesh(dat0[1],dat0[2],np.sum(diff,axis=0)/np.sum(dz))
plt.colorbar()

for i in range(0,xmax):
	for j in range(0,ymax):
		diff[:,j,i] = p0[8][:,j,i]*dz
plt.figure(5)
plt.title("vertically averaged vtheta")
plt.pcolormesh(dat0[1],dat0[2],np.sum(diff,axis=0)/np.sum(dz))
plt.colorbar()

for i in range(0,xmax):
	for k in range(0,zmax):
		diff[k,:,i] = p0[4][k,:,i]
plt.figure(6)
plt.title("density r-z slice")
plt.pcolormesh(dat0[1],np.pi/2.0-dat0[3],np.log10(0.5*(diff[:,ymax//2,:]+diff[:,ymax//2-1,:])))
wid = 2.5*0.035
plt.xlim(1.0-wid,1.0+wid)
plt.ylim(0.0,2.0*wid)
plt.colorbar()

for i in range(0,xmax):
	for k in range(0,zmax):
		diff[k,:,i] = p0[8][k,:,i]
plt.figure(7)
plt.title("vz r-z slice")
plt.pcolormesh(dat0[1],np.pi/2.0-dat0[3],0.5*(diff[:,ymax//2,:]+diff[:,ymax//2-1,:]))
wid = 2.5*0.035
plt.xlim(1.0-wid,1.0+wid)
plt.ylim(0.0,2.0*wid)
plt.colorbar()

plt.show()
