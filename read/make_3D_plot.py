import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd


label = 'h35_p14E_PPM'
xmax = 192+96
ymax = 672+96
zmax = 96+32

h = 0.035

frame = 11
dat0 = rd.load_3D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,zmax,label,0)
p0 = rd.load_3D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,zmax,label,frame)

dr = dat0[1][1:xmax+1]-dat0[1][0:xmax]
dp = dat0[2][1:ymax+1]-dat0[2][0:ymax]
dz = dat0[3][1:zmax+1]-dat0[3][0:zmax]

rc = 0.5*(dat0[1][1:xmax+1]+dat0[1][0:xmax])
pc = 0.5*(dat0[2][1:ymax+1]+dat0[2][0:ymax])
zc = 0.5*(dat0[3][1:zmax+1]+dat0[3][0:zmax])

ra = (dat0[1]-1.0)/h
pa = (dat0[2]-np.pi)/h
za = (np.pi/2.0-dat0[3])/h

plt.figure(0)
plt.plot((rc-1.0)/h,dr,"r")
plt.plot((pc-np.pi)/h,dp,"b")
plt.plot((zc-np.pi/2.0)/h,dz,"g")

diff = dat0[4]

for i in range(0,xmax):
	for j in range(0,ymax):
		diff[:,j,i] = p0[4][:,j,i]*dz
plt.figure(1)
plt.title("vertically integrated density")
plt.pcolormesh(ra,pa,np.log10(np.sum(diff,axis=0)/0.035/np.sqrt(np.pi/2.0)))
plt.colorbar()
wid = 2.0
plt.xlim(-wid,wid)
plt.ylim(-wid,wid)
plt.gca().set_aspect('equal')

for i in range(0,xmax):
	for j in range(0,ymax):
		diff[:,j,i] = p0[5][:,j,i]*dz
plt.figure(2)
plt.title("vertically integrated pressure")
plt.pcolormesh(ra,pa,np.log10(np.sum(diff,axis=0)/0.035/np.sqrt(np.pi/2.0)))
plt.colorbar()
wid = 2.0
plt.xlim(-wid,wid)
plt.ylim(-wid,wid)
plt.gca().set_aspect('equal')

for i in range(0,xmax):
	for j in range(0,ymax):
		diff[:,j,i] = p0[6][:,j,i]*p0[4][:,j,i]
plt.figure(3)
plt.title("vertically averaged vr")
plt.pcolormesh(ra,pa,np.sum(diff,axis=0)/np.sum(p0[4],axis=0))
plt.colorbar()
wid = 2.0
plt.xlim(-wid,wid)
plt.ylim(-wid,wid)
plt.gca().set_aspect('equal')

for i in range(0,xmax):
	for j in range(0,ymax):
		diff[:,j,i] = p0[7][:,j,i]*p0[4][:,j,i]
plt.figure(4)
plt.title("vertically averaged vphi")
plt.pcolormesh(ra,pa,np.sum(diff,axis=0)/np.sum(p0[4],axis=0))
plt.colorbar()
wid = 2.0
plt.xlim(-wid,wid)
plt.ylim(-wid,wid)
plt.gca().set_aspect('equal')

for i in range(0,xmax):
	for j in range(0,ymax):
		diff[:,j,i] = p0[8][:,j,i]*p0[4][:,j,i]
plt.figure(5)
plt.title("vertically averaged vtheta")
plt.pcolormesh(ra,pa,np.sum(diff,axis=0)/np.sum(p0[4],axis=0))
plt.colorbar()
wid = 2.0
plt.xlim(-wid,wid)
plt.ylim(-wid,wid)
plt.gca().set_aspect('equal')

for i in range(0,xmax):
	for k in range(0,zmax):
		diff[k,:,i] = p0[4][k,:,i]
plt.figure(60)
plt.title("density r-z slice")
plt.pcolormesh(ra,za,np.log10(0.5*(diff[:,ymax//2,:]+diff[:,ymax//2-1,:])))
wid = 2.0
plt.xlim(-wid,wid)
plt.ylim(0.0,2.0*wid)
plt.colorbar()
plt.gca().set_aspect('equal')

plt.figure(61)
plt.title("density r-z slice")
plt.pcolormesh(ra,za,np.log10(0.5*(diff[:,ymax//2,:]+diff[:,ymax//2-1,:])))
wid /= 2.0
plt.xlim(-wid,wid)
plt.ylim(0.0,2.0*wid)
plt.colorbar()
plt.gca().set_aspect('equal')

plt.figure(62)
plt.title("density r-z slice")
plt.pcolormesh(ra,za,np.log10(0.5*(diff[:,ymax//2,:]+diff[:,ymax//2-1,:])))
wid /= 2.0
plt.xlim(-wid,wid)
plt.ylim(0.0,2.0*wid)
plt.colorbar()
plt.gca().set_aspect('equal')

plt.figure(63)
plt.title("density r-z slice")
plt.pcolormesh(ra,za,np.log10(0.5*(diff[:,ymax//2,:]+diff[:,ymax//2-1,:])))
wid /= 2.0
plt.xlim(-wid,wid)
plt.ylim(0.0,2.0*wid)
plt.colorbar()
plt.gca().set_aspect('equal')

for i in range(0,xmax):
	for k in range(0,zmax):
		diff[k,:,i] = p0[8][k,:,i]/h
plt.figure(7)
plt.title("vz r-z slice")
plt.pcolormesh(ra,za,0.5*(diff[:,ymax//2,:]+diff[:,ymax//2-1,:]),cmap="seismic")
wid = 2.0*h
plt.xlim(-2.0,2.0)
plt.ylim(0.0,4.0)
plt.clim(-10.0,10.0)
plt.colorbar()
plt.gca().set_aspect('equal')


plt.show()
