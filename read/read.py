import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
from matplotlib import cm
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import math
import pickle
import os.path as pat

twopi = 2.0*np.pi

def check_file_exist(path, imax, jmax, label, num):
	fname = path+'binary_'
	fname = fname + str(imax)+'x'+str(jmax)+'_'
	fname = fname + label + '_'
	fname = fname + frame_num(num)
	return pat.isfile(fname)

def cell_center(xa):
	N = len(xa)-1
	xc = np.arange(N,dtype=np.float64)

	for i in range(N):
		xc[i] = (xa[i+1]+xa[i])/2.0
	return xc

def load_data(fname):
	f = open(fname+'.pckl','rb')
	out = pickle.load(f)
	f.close()
	return out

def save_data(dat, fname):
	f = open(fname+'.pckl','wb')
	pickle.dump(dat,f)
	f.close()
	return fname+'.pckl'

def frame_num(i):
	snum = str(i)
	if len(snum)==1:   snum = '0000'+snum
	elif len(snum)==2: snum = '000'+snum
	elif len(snum)==3: snum = '00'+snum
	elif len(snum)==4: snum = '0'+snum
	return snum

def periodic_bound(x,xmin,xmax):
	while x<xmin: x += xmax-xmin
	while x>=xmax: x -= xmax-xmin
	return x

def ceil_2pow(val):
	tar = 1
	while tar<val: tar *= 2
	return tar

def bin_locate(val,arr):
	L = len(arr)-1
	if val<arr[0] or val>arr[L]: return None
	if val==arr[L]: return L

	N = ceil_2pow(L)//2
	mid = N
	while N>1:
		N = N//2
		if mid >= L:       mid -= N 
		elif arr[mid]>val: mid -= N
		else:              mid += N
	if arr[mid]>val: mid -= 1
	return mid

def interpolate_2D_periodicY(p,x,y,z,ymin,ymax):
	px = p[0]
	py = periodic_bound(p[1],ymin,ymax)

	jmax = len(y)
	ic = bin_locate(px,x)
	jc = bin_locate(py,y)
	if ic==None: return None
	
	if ic==len(x)-1: 
		A = 0.0
		ic1 = ic
	else:
		A = (px-x[ic])/(x[ic+1]-x[ic])
		ic1 = ic+1

	if jc==jmax-1 or jc==None:
		jc = len(y)-1
		jc1 = 0
		if   py<y[0]:       B = 0.5 + 0.5*(py-ymin)/(y[0]-ymin)
		elif py>=y[jmax-1]: B = 0.5*(py-y[jmax-1])/(ymax-y[jmax-1])
		else: return None
	else:
		B = (py-y[jc])/(y[jc+1]-y[jc])
		jc1 = jc+1

	Q1 = A*B
	Q2 = (1.0-A)*B
	Q3 = (1.0-A)*(1.0-B)
	Q4 = A*(1.0-B)

	val = Q1*z[jc1,ic1] + Q2*z[jc1,ic] + Q3*z[jc,ic] + Q4*z[jc,ic1]
	return val

def interpolate_2D(p,x,y,z):
	ic = bin_locate(p[0],x)
	jc = bin_locate(p[1],y)
	if ic==None or jc==None: return None
	
	if ic==len(x)-1: 
		A = 0.0
		ic1 = ic
	else:
		A = (p[0]-x[ic])/(x[ic+1]-x[ic])
		ic1 = ic+1

	if jc==len(y)-1:
		B = 0.0
		jc1 = jc
	else:	
		B = (p[1]-y[jc])/(y[jc+1]-y[jc])
		jc1 = jc+1

	Q1 = A*B
	Q2 = (1.0-A)*B
	Q3 = (1.0-A)*(1.0-B)
	Q4 = A*(1.0-B)

	val = Q1*z[jc1,ic1] + Q2*z[jc1,ic] + Q3*z[jc,ic] + Q4*z[jc,ic1]
	return val

def get_vorticity_z(xc,yc,vx,vy):
	imax = len(xc)
	jmax = len(yc)

	vor = np.arange(imax*jmax, dtype=np.float64).reshape(jmax,imax)

	for j in range(jmax):
		for i in range(imax):
			vor[j,i] = 0.0

	for j in range(1,jmax-1):
		for i in range(1,imax-1):
			vor[j,i] = (vy[j,i+1]-vy[j,i-1])/(xc[i+1]-xc[i-1])
			vor[j,i]-= (vx[j+1,i]-vx[j-1,i])/(yc[j+1]-yc[j-1])
	return vor

def load_1D_data(path, imax, label, num):
	fname = path+'binary_'
	fname = fname + str(imax)+'_'
	fname = fname + label + '_'
	fname = fname + frame_num(num)
	
	dat = np.fromfile(fname, dtype=np.float64)

	start = 0
	end = 1
	t = dat[0]

	start = end
	end = start+imax+1
	xa = dat[start:end]

	xc = np.arange(imax, dtype=np.float64)
	for i in range(0,imax): xc[i] = 0.5*(xa[i+1]+xa[i])

	start = end
	end = start+imax
	ro = dat[start:end]

	start = end
	end = start+imax
	pr = dat[start:end]

	start = end
	end = start+imax
	vx = dat[start:end]

	return (t,xc,ro,pr,vx)

def load_1D_dust_data(path, imax, label, num):
	fname = path+'binary_'
	fname = fname + str(imax)+'x1x1_'
	fname = fname + label + '_'
	fname = fname + frame_num(num)
	
	dat = np.fromfile(fname, dtype=np.float64)

	start = 0
	end = 1
	t = dat[0]

	start = end
	end = start+imax+1
	xa = dat[start:end]

	xc = np.arange(imax, dtype=np.float64)
	for i in range(0,imax): xc[i] = 0.5*(xa[i+1]+xa[i])

	start = end
	end = start+imax
	ro = dat[start:end]

	start = end
	end = start+imax
	pr = dat[start:end]

	start = end
	end = start+imax
	vx = dat[start:end]

	start = end
	end = start+imax
	ro_d = dat[start:end]

	start = end
	end = start+imax
	vx_d = dat[start:end]

	return (t,xc,ro,pr,vx,ro_d,vx_d)

def load_2D_data(path, imax, jmax, label, num):
	fname = path+'binary_'
	fname = fname + str(imax)+'x'+str(jmax)+'_'
	fname = fname + label + '_'
	fname = fname + frame_num(num)
	
	dat = np.fromfile(fname, dtype=np.float64)

	start = 0
	end = 1
	time = dat[0]

	start = end
	end = start+imax+1
	x = dat[start:end]

	start = end
	end = start+jmax+1
	y = dat[start:end]

	start = end
	end = start+imax*jmax
	ro = np.reshape(dat[start:end], (-1,imax))

	start = end
	end = start+imax*jmax
	pr = np.reshape(dat[start:end], (-1,imax))

	start = end
	end = start+imax*jmax
	vx = np.reshape(dat[start:end], (-1,imax))

	start = end
	end = start+imax*jmax
	vy = np.reshape(dat[start:end], (-1,imax))

	return (time,x,y,ro,pr,vx,vy)

def load_2D_dust_data(path, imax, jmax, label, num):
	fname = path+'binary_'
	fname = fname + str(imax)+'x'+str(jmax)+'_'
	fname = fname + label + '_'
	fname = fname + frame_num(num)
	
	dat = np.fromfile(fname, dtype=np.float64)

	start = 0
	end = 1
	time = dat[0]

	start = end
	end = start+imax+1
	x = dat[start:end]

	start = end
	end = start+jmax+1
	y = dat[start:end]

	start = end
	end = start+imax*jmax
	ro = np.reshape(dat[start:end], (-1,imax))

	start = end
	end = start+imax*jmax
	pr = np.reshape(dat[start:end], (-1,imax))

	start = end
	end = start+imax*jmax
	vx = np.reshape(dat[start:end], (-1,imax))

	start = end
	end = start+imax*jmax
	vy = np.reshape(dat[start:end], (-1,imax))

	start = end
	end = start+imax*jmax
	ro_d = np.reshape(dat[start:end], (-1,imax))

	start = end
	end = start+imax*jmax
	vx_d = np.reshape(dat[start:end], (-1,imax))

	start = end
	end = start+imax*jmax
	vy_d = np.reshape(dat[start:end], (-1,imax))

	return (time,x,y,ro,pr,vx,vy,ro_d,vx_d,vy_d)
	
def load_3D_data(path, imax, jmax, kmax, label, num):
	fname = path+'binary_'
	fname = fname + str(imax)+'x'+str(jmax)+'x'+str(kmax)+'_'
	fname = fname + label + '_'
	fname = fname + frame_num(num)
	
	dat = np.fromfile(fname, dtype=np.float64)

	start = 0
	end = 1
	time = dat[0]

	start = end
	end = start+imax+1
	x = dat[start:end]

	start = end
	end = start+jmax+1
	y = dat[start:end]
	
	start = end
	end = start+kmax+1
	z = dat[start:end]

	start = end
	end = start+imax*jmax*kmax
	ro = np.reshape(dat[start:end], (-1,jmax,imax))

	start = end
	end = start+imax*jmax*kmax
	pr = np.reshape(dat[start:end], (-1,jmax,imax))

	start = end
	end = start+imax*jmax*kmax
	vx = np.reshape(dat[start:end], (-1,jmax,imax))

	start = end
	end = start+imax*jmax*kmax
	vy = np.reshape(dat[start:end], (-1,jmax,imax))
	
	start = end
	end = start+imax*jmax*kmax
	vz = np.reshape(dat[start:end], (-1,jmax,imax))

	return (time,x,y,z,ro,pr,vx,vy,vz)

def load_3D_dust_data(path, imax, jmax, kmax, label, num):
	fname = path+'binary_'
	fname = fname + str(imax)+'x'+str(jmax)+'x'+str(kmax)+'_'
	fname = fname + label + '_'
	fname = fname + frame_num(num)
	
	dat = np.fromfile(fname, dtype=np.float64)

	start = 0
	end = 1
	time = dat[0]

	start = end
	end = start+imax+1
	x = dat[start:end]

	start = end
	end = start+jmax+1
	y = dat[start:end]
	
	start = end
	end = start+kmax+1
	z = dat[start:end]

	start = end
	end = start+imax*jmax*kmax
	ro = np.reshape(dat[start:end], (-1,jmax,imax))

	start = end
	end = start+imax*jmax*kmax
	pr = np.reshape(dat[start:end], (-1,jmax,imax))

	start = end
	end = start+imax*jmax*kmax
	vx = np.reshape(dat[start:end], (-1,jmax,imax))

	start = end
	end = start+imax*jmax*kmax
	vy = np.reshape(dat[start:end], (-1,jmax,imax))
	
	start = end
	end = start+imax*jmax*kmax
	vz = np.reshape(dat[start:end], (-1,jmax,imax))

	start = end
	end = start+imax*jmax*kmax
	ro_d = np.reshape(dat[start:end], (-1,jmax,imax))

	start = end
	end = start+imax*jmax*kmax
	vx_d = np.reshape(dat[start:end], (-1,jmax,imax))

	start = end
	end = start+imax*jmax*kmax
	vy_d = np.reshape(dat[start:end], (-1,jmax,imax))
	
	start = end
	end = start+imax*jmax*kmax
	vz_d = np.reshape(dat[start:end], (-1,jmax,imax))

	return (time,x,y,z,ro,pr,vx,vy,vz,ro_d,vx_d,vy_d,vz_d)

def surface_density_half_disk(i,j,z,ro):
	dum = 0.0
	for k in range(0,len(z)-1):
		dum += ro[k,j,i]*(z[k+1]-z[k])
	return dum/(0.05*np.sqrt(np.pi/2.0))

def make_cart_plot(time,x,y,z,val):
	fig,ax = plt.subplots(figsize=(7,7))
	
	pres = 512
	imax = len(x)-1
	jmax = len(y)-1
	kmax = len(z)-1
	
	rad = np.zeros(imax)
	azi = np.zeros(jmax)
	sig2D = np.reshape(np.zeros(imax*jmax), (-1,imax))

	box_size = 7.0
	cx = np.zeros(pres+1)
	cy = np.zeros(pres+1)
	dd = 0.5*box_size/pres
	aa = np.reshape(np.zeros(pres*pres), (-1,pres))
	
	for i in range(0,pres+1):
		cx[i] = (box_size*i)/pres - 0.5*box_size;
		cy[i] = (box_size*i)/pres - 0.5*box_size;

	for i in range(0,imax):
		rad[i] = (x[i+1]+x[i])/2.0
		for j in range(0,jmax):
			azi[j] = (y[j+1]+y[j])/2.0
			sig2D[j,i] = surface_density_half_disk(i,j,z,val)*rad[i]

	for i in range(0,pres):
		for j in range(0,pres):
			dumx = np.sqrt((cx[i]+dd)*(cx[i]+dd) + (cy[j]+dd)*(cy[j]+dd))
			dumy = np.arctan((cy[j]+dd)/(cx[i]+dd))
			if cx[i]<0.0: dumy += np.pi
			elif cy[j]<0.0: dumy += 2.0*np.pi
			aa[j,i] = interpolate_2D_periodicY((dumx,dumy),rad,azi,sig2D,0.0,2.0*np.pi)
		
	mesh1 = ax.pcolormesh(cx,cy,np.log10(aa),cmap='magma',vmin=-0.6, vmax=0.4)
	ax.set_aspect('equal')
	
	cbar = fig.colorbar(mesh1, ax=ax, shrink=0.65)
	cbar.ax.tick_params(labelsize=18, width=2)
	
	stime = str(round(time/(2.0*np.pi),2))
	ax.set_title('t = '+stime, fontsize=20)
	ax.set_xlabel('x', fontsize=18)
	ax.set_ylabel('y', fontsize=18)
	ax.set_xticks(np.arange(-3.0, 3.1, step=1.0))
	ax.set_yticks(np.arange(-3.0, 3.1, step=1.0))
	ax.set_xlim(-0.5*box_size,0.5*box_size)
	ax.set_ylim(-0.5*box_size,0.5*box_size)
	ax.tick_params(axis='both', labelsize=18, width=2)
	fig.tight_layout()

	return fig,ax
	
def make_polar_plot_2D(fig,ax,time,x,y,val,vmin,vmax):
	imax = len(x)-1
	jmax = len(y)-1

	rad = np.zeros(imax)
	for i in range(0,imax):
		rad[i] = (x[i+1]+x[i])/2.0

	if val.min()>0.0:
		mesh1 = ax.pcolormesh(x,y,np.log10(val),cmap='magma',vmin=vmin, vmax=vmax)
	else:
		vrange = max(val.max(),-val.min())
		mesh1 = ax.pcolormesh(x,y,val,cmap='PuOr',vmin=-vrange,vmax=vrange)

	cbar = fig.colorbar(mesh1, ax=ax, shrink=0.65)
	cbar.ax.tick_params(labelsize=18, width=2)
	
	stime = str(round(time/(2.0*np.pi),2))
	ax.set_title('t = '+stime, fontsize=20)
	ax.set_xlabel('r', fontsize=18)
	ax.set_ylabel(r'$\phi$', fontsize=18)
	ax.tick_params(axis='both', labelsize=18, width=2)

	return mesh1

def make_cart_plot_2D(fig,ax,time,x,y,val,vmin,vmax):
	pres = 512
	imax = len(x)-1
	jmax = len(y)-1
	
	rad = np.zeros(imax)
	azi = np.zeros(jmax)

	box_size = 6.0
	cx = np.zeros(pres+1)
	cy = np.zeros(pres+1)
	dd = 0.5*box_size/pres
	aa = np.reshape(np.zeros(pres*pres), (-1,pres))
	
	for i in range(0,pres+1):
		cx[i] = (box_size*i)/pres - 0.5*box_size;
		cy[i] = (box_size*i)/pres - 0.5*box_size;

	for i in range(0,imax):
		rad[i] = (x[i+1]+x[i])/2.0
		for j in range(0,jmax):
			azi[j] = (y[j+1]+y[j])/2.0

	for i in range(0,pres):
		for j in range(0,pres):
			dumx = np.sqrt((cx[i]+dd)*(cx[i]+dd) + (cy[j]+dd)*(cy[j]+dd))
			dumy = np.arctan((cy[j]+dd)/(cx[i]+dd))
			if cx[i]<0.0: dumy += np.pi
			elif cy[j]<0.0: dumy += 2.0*np.pi
			if dumx>=rad.max(): aa[j,i]=np.power(10.0,vmin)
			else: aa[j,i] = interpolate_2D_periodicY((dumx,dumy),rad,azi,val,0.0,2.0*np.pi)
			

	if val.min()>0.0:
		mesh1 = ax.pcolormesh(cx,cy,np.log10(aa),cmap='magma',vmin=vmin, vmax=vmax)
	else:
		vrange = max(val.max(),-val.min())
		mesh1 = ax.pcolormesh(cx,cy,aa,cmap='PuOr',vmin=-vrange,vmax=vrange)

	ax.set_aspect('equal')
	ax.set_xlabel('x', fontsize=12)
	ax.set_ylabel('y', fontsize=12)
	ax.set_xticks(np.arange(-5.0, 5.1, step=1.0))
	ax.set_yticks(np.arange(-5.0, 5.1, step=1.0))
	ax.set_xlim(-0.5*box_size,0.5*box_size)
	ax.set_ylim(-0.5*box_size,0.5*box_size)
	ax.tick_params(axis='both', labelsize=12, width=2)

	return mesh1

def make_polar_plot_3D(fig,ax,x,y,z,val,vmin,vmax):
	imax = len(x)-1
	jmax = len(y)-1
	kmax = len(z)-1
	
	rad = np.zeros(imax)
	azi = np.zeros(jmax)
	sig2D = np.reshape(np.zeros(imax*jmax), (-1,imax))

	for i in range(0,imax):
		rad[i] = (x[i+1]+x[i])/2.0
		for j in range(0,jmax):
			azi[j] = (y[j+1]+y[j])/2.0
			sig2D[j,i] = surface_density_half_disk(i,j,z,val)*rad[i]

	if sig2D.min()>0.0:
		mesh1 = ax.pcolormesh(x,y,np.log10(sig2D),cmap='magma',vmin=vmin, vmax=vmax)
	else:
		vrange = max(sig2D.max(),-sig2D.min())
		mesh1 = ax.pcolormesh(x,y,sig2D,cmap='PuOr',vmin=-vrange,vmax=vrange)

	ax.set_xlim(0.4,4.0)
	ax.set_aspect(0.67*(x.max()-x.min())/(y.max()-y.min()))
	ax.set_xlabel('r', fontsize=16)
	ax.set_ylabel('$\phi$', fontsize=16)
	ax.tick_params(axis='both', labelsize=16)

	return mesh1

def make_rz_plot_3D(fig,ax,x,y,z,val,vmin,vmax):
	imax = len(x)-1
	jmax = len(y)-1
	kmax = len(z)-1
	
	rad = np.zeros(imax)
	pol = np.zeros(kmax)
	zp  = np.zeros(kmax+1)
	den2D = np.reshape(np.zeros(imax*kmax), (-1,imax))

	for i in range(0,imax):
		rad[i] = (x[i+1]+x[i])/2.0
		for k in range(0,kmax):
			pol[k] = (z[k+1]+z[k])/2.0
			zp [k] = 0.5*np.pi-z[k];
			den2D[k,i] = val[k,:,i].mean()
	zp [kmax] = 0.5*np.pi-z[kmax];

	mesh1 = ax.pcolormesh(x,zp,np.log10(den2D),cmap='magma',vmin=vmin, vmax=vmax)

	ax.set_xlim(0.4,4.0)
	ax.set_aspect(0.67*(x.max()-x.min())/(zp.max()-zp.min()))
	ax.set_xlabel('r', fontsize=16)
	ax.set_ylabel('$\pi/2 - \\theta$', fontsize=16)
	ax.tick_params(axis='both', labelsize=16)

	return mesh1

def make_figure_2D(dat):
	fig, axs = plt.subplots(1, 2, figsize=(11,5))
	m1=make_cart_plot_2D(fig,axs[0],dat[0],dat[1],dat[2],dat[3],-1.5,0.5)#,-3.0,1.0)
	m2=make_cart_plot_2D(fig,axs[1],dat[0],dat[1],dat[2],dat[7]/100.0,-4.0,0.5)#,-7.0,1.0)
	plt.figtext(0.41,0.88,'t = '+str(round(dat[0]/(2.0*np.pi),2)/1000.0)+' kyr', fontsize=18)
	plt.figtext(0.24,0.84,'Gas', fontsize=16)
	plt.figtext(0.67,0.84,'Dust', fontsize=16)
	cbar1 = fig.colorbar(m1, ax=axs[0], shrink=0.6, pad=0.02)
	cbar1.ax.tick_params(labelsize=12, width=1)

	cbar2 = fig.colorbar(m2, ax=axs[1], shrink=0.6, pad=0.02)
	cbar2.ax.tick_params(labelsize=12, width=1)
	return fig,axs

def make_figure_2D_nd(dat):
	fig, axs = plt.subplots(1, 2, figsize=(11,5))
	m1=make_cart_plot_2D(fig,axs[0],dat[0],dat[1],dat[2],dat[3],-0.75,0.5)
	m2=make_cart_plot_2D(fig,axs[1],dat[0],dat[1],dat[2],dat[4],-4.25,-2.75)
	plt.figtext(0.41,0.88,'t = '+str(round(dat[0]/(2.0*np.pi),2))+' yr', fontsize=18)
	plt.figtext(0.24,0.84,'density', fontsize=16)
	plt.figtext(0.67,0.84,'pressure', fontsize=16)
	cbar1 = fig.colorbar(m1, ax=axs[0], shrink=0.6, pad=0.02)
	cbar1.ax.tick_params(labelsize=12, width=1)

	cbar2 = fig.colorbar(m2, ax=axs[1], shrink=0.6, pad=0.02)
	cbar2.ax.tick_params(labelsize=12, width=1)
	return fig,axs

def make_figure_2D_nd(dat,dat0):
	fig, axs = plt.subplots(1, 2, figsize=(11,5))
	m1=make_cart_plot_2D(fig,axs[0],dat[0],dat[1],dat[2],dat[3]-dat0[3],-0.75,0.5)
	m2=make_cart_plot_2D(fig,axs[1],dat[0],dat[1],dat[2],dat[4]-dat0[4],-4.25,-2.75)
	plt.figtext(0.41,0.88,'t = '+str(round(dat[0]/(2.0*np.pi),2))+' yr', fontsize=18)
	plt.figtext(0.24,0.84,'density', fontsize=16)
	plt.figtext(0.67,0.84,'pressure', fontsize=16)
	cbar1 = fig.colorbar(m1, ax=axs[0], shrink=0.6, pad=0.02)
	cbar1.ax.tick_params(labelsize=12, width=1)

	cbar2 = fig.colorbar(m2, ax=axs[1], shrink=0.6, pad=0.02)
	cbar2.ax.tick_params(labelsize=12, width=1)
	return fig,axs

def make_figure_3D(dat):
	fig, axs = plt.subplots(2, 2, figsize=(13,8), sharex=True)
	m1=make_polar_plot_3D(fig,axs[0,0],dat[1],dat[2],dat[3],dat[4],-1.5,0.5)
	m3=make_polar_plot_3D(fig,axs[1,0],dat[1],dat[2],dat[3],dat[9]/100.0,-4.0,0.5)
	m2=make_rz_plot_3D(fig,axs[0,1],dat[1],dat[2],dat[3],dat[4],-1.5,0.5)
	m4=make_rz_plot_3D(fig,axs[1,1],dat[1],dat[2],dat[3],dat[9]/100.0,-4.0,0.5)
	#plt.figtext(0.43,0.94,'t = '+str(round(dat[0]/(2.0*np.pi),2))+' orbits', fontsize=18)
	#plt.figtext(0.26,0.9,'Gas', fontsize=16)
	#plt.figtext(0.69,0.9,'Dust', fontsize=16)

	axs[0,0].set_xlabel('')
	axs[0,1].set_xlabel('')
	#axs[0,1].set_ylabel('')
	#axs[1,1].set_ylabel('')

	fig.text(0.05,0.73,'Gas', fontsize=16, rotation=90)
	fig.text(0.05,0.26,'Dust', fontsize=16, rotation=90)

	fig.tight_layout()

	cbar1 = fig.colorbar(m2, ax=axs[0,:].ravel().tolist(), shrink=0.6, location='right', pad=0.02)
	cbar1.ax.tick_params(labelsize=14, width=1)

	cbar2 = fig.colorbar(m4, ax=axs[1,:].ravel().tolist(), shrink=0.6, location='right', pad=0.02)
	cbar2.ax.tick_params(labelsize=14, width=1)

	fig.savefig('../images/gas_dust_3D.png', dpi=300, bbox_inches='tight')
	return fig,axs

def make_dM_dr_plot(x,y,z,tr0,tr1,tr2,tr3,ratio):
	fig,ax = plt.subplots(figsize=(8,6))
	
	imax = tr0.shape[2]
	rad = np.zeros(imax)
	sig0 = np.zeros(imax)
	sig1 = np.zeros(imax)
	sig2 = np.zeros(imax)
	sig3 = np.zeros(imax)
	
	for i in range(0,imax):
		rad[i] = (x[i+1]+x[i])/2.0
		for k in range(0,kmax):
			dz = z[k+1]-z[k]
			sig0[i] += (tr0[k,:,i]*dz*rad[i]).mean()/(0.05*np.sqrt(np.pi/2.0))
			sig1[i] += (tr1[k,:,i]*dz*rad[i]).mean()/(0.05*np.sqrt(np.pi/2.0))
			sig2[i] += (tr2[k,:,i]*dz*rad[i]).mean()/(0.05*np.sqrt(np.pi/2.0))
			sig3[i] += (tr3[k,:,i]*dz*rad[i]).mean()/(0.05*np.sqrt(np.pi/2.0))
		sig0[i] *= rad[i]
		sig1[i] *= rad[i]
		sig2[i] *= rad[i]
		sig3[i] *= rad[i]

	ax.plot(rad,sig0,'--',color='grey')
	ax.plot(rad,sig1,'r')
	ax.plot(rad,sig2,'g')
	ax.plot(rad,sig3,'b')
	ax.set_xlim(1.0/ratio,2.0*ratio)
	ax.set_ylim(0.4,1.4)
	
	return fig,ax
	
def make_M_r_plot(imax,kmax,x,y,z,tr0,tr1,tr2,tr3,ratio):
	fig,ax = plt.subplots(figsize=(8,6))
	
	imax = tr0.shape[2]
	rad = np.zeros(imax)
	sig0 = np.zeros(imax)
	sig1 = np.zeros(imax)
	sig2 = np.zeros(imax)
	sig3 = np.zeros(imax)
	
	for i in range(0,imax):
		rad[i] = (x[i+1]+x[i])/2.0
		for k in range(0,kmax):
			dz = z[k+1]-z[k]
			sig0[i] += (tr0[k,:,i]*dz*rad[i]).mean()/(0.05*np.sqrt(np.pi/2.0))
			sig1[i] += (tr1[k,:,i]*dz*rad[i]).mean()/(0.05*np.sqrt(np.pi/2.0))
			sig2[i] += (tr2[k,:,i]*dz*rad[i]).mean()/(0.05*np.sqrt(np.pi/2.0))
			sig3[i] += (tr3[k,:,i]*dz*rad[i]).mean()/(0.05*np.sqrt(np.pi/2.0))
		sig0[i] *= rad[i]
		sig1[i] *= rad[i]
		sig2[i] *= rad[i]
		sig3[i] *= rad[i]
		
	for i in range(0,imax):
		dr = x[i+1]-x[i]
		sig0[i] *= dr
		sig1[i] *= dr
		sig2[i] *= dr
		sig3[i] *= dr

	istart = bin_locate(1.0/ratio,rad)
	#iend = bin_locate(2.0*ratio,rad)

	for i in range(istart,imax):
		sig0[i] += sig0[i-1]
		sig1[i] += sig1[i-1]
		sig2[i] += sig2[i-1]
		sig3[i] += sig3[i-1]
	sig0[istart-1] = sig0[istart]/1000.0
	sig1[istart-1] = sig0[istart]/1000.0
	sig2[istart-1] = sig0[istart]/1000.0
	sig3[istart-1] = sig0[istart]/1000.0
	
	print(sig3[istart-1])
	print(sig3[istart])
	print(sig0[istart-1])
	print(sig0[istart])

	ax.plot(x[istart:],(sig0[istart-1:]/sig0[istart-1:])-1.0,'--',color='grey')
	ax.plot(x[istart:],(sig1[istart-1:]/sig0[istart-1:])-1.0,'r')
	ax.plot(x[istart:],(sig2[istart-1:]/sig0[istart-1:])-1.0,'g')
	ax.plot(x[istart:],(sig3[istart-1:]/sig0[istart-1:])-1.0,'b')
	ax.set_xlim(1.0/ratio,2.0*ratio)
	#ax.set_ylim(-0.06,0.04)
	
	return fig,ax

def make_PPD_plots():
	imax = 480
	jmax = 1280
	kmax = 40
	path = '/data/fung/'
	label = 'c50_p42_s200_r100_i45_o500_ad_rot_fargo_bc10_PEM'

	for i in range(0,71):
		(t,x,y,z,ro,pr,vx,vy,vz,tr)=load_3D_data(path,imax,jmax,kmax,label,i)
		f1,a1=make_cart_plot(t,x,y,z,ro)
		f1.savefig('../images/frame_'+frame_num(i)+'_den.png', dpi=300, bbox_inches='tight')
		plt.close()
	return

def make_gas_dust_movie_2D(istart,iend):
	for i in range(istart,iend+1):
		dat = load_2D_dust_data('/data/fung/',304,754,'c50_p21_s5000_r100_i40_o500_b10_rot_fargo_PEMEU', i)
		f1,axes=make_figure_2D(dat)
		f1.savefig('../images/2D_iso_'+frame_num(i)+'.png', dpi=300, bbox_inches='tight')
		plt.close('all')
		print('frame',i,'done.')
	return

def make_gas_dust_movie_3D(istart,iend):
	for i in range(istart,iend+1):
		dat = load_3D_dust_data('/data/fung/',304,754,36,'c50_p21_s250_r100_i40_o500_ad_rot_fargo_bc10_dr_PEMEU',i)
		f1,axes=make_figure_3D(dat)
		f1.savefig('../images/3D_'+frame_num(i)+'.png', dpi=300, bbox_inches='tight')
		plt.close('all')
		print('frame',i,'done.')
	return

def ssa_plot():
	f1,a1=plt.subplots()

	r1_ssa = load_2D_dust_data('/data/fung/',304,754,'c50_i40_o500_b10_fargo_PEMEU', 1)
	r2_ssa = load_2D_dust_data('/data/fung/',304*2,754*2,'c50_i40_o500_b10_fargo_PEMEU', 1)
	r3_ssa = load_2D_dust_data('/data/fung/',304*4,754*4,'c50_i40_o500_b10_fargo_PEMEU', 1)
	r4_ssa = load_2D_dust_data('/data/fung/',304*8,754*8,'c50_i40_o500_b10_fargo_PEMEU', 1)

	r1_sa = load_2D_dust_data('/data/fung/',304,754,'c50_i40_o500_b10_fargo_tmp_PEMEU', 1)
	r2_sa = load_2D_dust_data('/data/fung/',304*2,754*2,'c50_i40_o500_b10_fargo_tmp_PEMEU', 1)
	r3_sa = load_2D_dust_data('/data/fung/',304*4,754*4,'c50_i40_o500_b10_fargo_tmp_PEMEU', 1)

	std = np.zeros(304)
	rad = np.zeros(304)
	for i in range(0,304):
		rad[i] = 0.5*(r1_sa[1][i]+r1_sa[1][i+1])
		std[i] = r4_ssa[8][:,i*8:i*8+8].mean()
	a1.plot(rad,np.abs(r1_sa[8][0,:]/std-1.0),'r:', label='SA; 304x754')

	std = np.zeros(304*2)
	rad = np.zeros(304*2)
	for i in range(0,304*2):
		rad[i] = 0.5*(r2_sa[1][i]+r2_sa[1][i+1])
		std[i] = r4_ssa[8][:,i*4:i*4+4].mean()
	a1.plot(rad,np.abs(r2_sa[8][0,:]/std-1.0),'b:', label='SA; 608x1508')

	std = np.zeros(304*4)
	rad = np.zeros(304*4)
	for i in range(0,304*4):
		rad[i] = 0.5*(r3_ssa[1][i]+r3_ssa[1][i+1])
		std[i] = r4_ssa[8][:,i*2:i*2+2].mean()
	a1.plot(rad,np.abs(r3_sa[8][0,:]/std-1.0),'k:', label='SA; 1216x3016')

	std = np.zeros(304)
	rad = np.zeros(304)
	for i in range(0,304):
		rad[i] = 0.5*(r1_ssa[1][i]+r1_ssa[1][i+1])
		std[i] = r4_ssa[8][:,i*8:i*8+8].mean()
	a1.plot(rad,np.abs(r1_ssa[8][0,:]/std-1.0),'r', label='SSA; 304x754')

	std = np.zeros(304*2)
	rad = np.zeros(304*2)
	for i in range(0,304*2):
		rad[i] = 0.5*(r2_ssa[1][i]+r2_ssa[1][i+1])
		std[i] = r4_ssa[8][:,i*4:i*4+4].mean()
	a1.plot(rad,np.abs(r2_ssa[8][0,:]/std-1.0),'b', label='SSA; 608x1508')

	std = np.zeros(304*4)
	rad = np.zeros(304*4)
	for i in range(0,304*4):
		rad[i] = 0.5*(r3_ssa[1][i]+r3_ssa[1][i+1])
		std[i] = r4_ssa[8][:,i*2:i*2+2].mean()
	a1.plot(rad,np.abs(r3_ssa[8][0,:]/std-1.0),'k', label='SSA; 1216x3016')

	a1.set_yscale('log')
	a1.set_xlim(0.4,5.0)
	a1.set_ylim(0.3e-7,5.0e-2)
	a1.tick_params(axis='x', labelsize=16)
	a1.tick_params(axis='y', labelsize=16)
	a1.set_xlabel('radius [au]', fontsize=16)
	a1.set_ylabel('relative error in drift speed', fontsize=16)
	a1.legend(fontsize=12)

	f1.tight_layout()
	f1.savefig('../images/SSA_err.png', dpi=300, bbox_inches='tight')
	return f1

def gas_ring_2D():
	time = np.zeros(1001)
	mass = np.zeros(1001)
	dat = load_2D_dust_data('/data/fung/',304,754,'c50_p21_s5000_r100_i40_o500_ad_rot_fargo_bc10_dr_PEMEU', 0)
	rad = np.zeros(304)
	for i in range(0,304):
		rad[i] = 0.5*(dat[1][i]+dat[1][i+1])

	for n in range(0,1001):
		dat = load_2D_dust_data('/data/fung/',304,754,'c50_p21_s5000_r100_i40_o500_ad_rot_fargo_bc10_dr_PEMEU', n)
		time[n] = dat[0]
		mass[n] = 0.0
		for i in range(0,304):
			if rad[i]>1.0 and rad[i]<2.0:
				mass[n] += dat[3][:,i].mean()*rad[i]*(dat[1][i+1]-dat[1][i])
		print('At t=',time[n],' ring mass is',mass[n])
	return time, mass

def gas_ring_plot(f1, a1, time, ring):
	tmp = time/np.pi/2000.0
	a1.plot(tmp,ring,'k',label='ring mass')
	fit = np.zeros(tmp.shape)
	fit = 0.895*np.exp(-tmp/500.0)
	a1.plot(tmp,fit,'--',color='grey', label=r'$\propto$exp(-t/500kyr)')

	a1.legend(fontsize=16)
	a1.set_xlabel('time [kyr]', fontsize=16)
	a1.set_ylabel('normalized ring mass', fontsize=16)
	a1.tick_params(axis='x', labelsize=16)
	a1.tick_params(axis='y', labelsize=16)
	return

def ring_plots():
	f1,ax=plt.subplots(1, 2, figsize=(11,5))

	dat = load_2D_dust_data('/data/fung/',304,754,'c50_p21_s5000_r100_i40_o500_ad_rot_fargo_bc10_dr_PEMEU', 0)
	rad = np.zeros(304)
	val = np.zeros(304)
	for i in range(0,304):
		rad[i] = 0.5*(dat[1][i]+dat[1][i+1])
		val[i] = dat[3][:,i].mean()
	ax[0].plot(rad,val,'--',color='grey',label='t=0kyr')

	dat = load_2D_dust_data('/data/fung/',304,754,'c50_p21_s5000_r100_i40_o500_ad_rot_fargo_bc10_dr_PEMEU', 10)
	rad = np.zeros(304)
	val = np.zeros(304)
	for i in range(0,304):
		rad[i] = 0.5*(dat[1][i]+dat[1][i+1])
		val[i] = dat[3][:,i].mean()
	ax[0].plot(rad,val,'k',label='t=1kyr')

	dat = load_2D_dust_data('/data/fung/',304,754,'c50_p21_s5000_r100_i40_o500_ad_rot_fargo_bc10_dr_PEMEU', 100)
	rad = np.zeros(304)
	val = np.zeros(304)
	for i in range(0,304):
		rad[i] = 0.5*(dat[1][i]+dat[1][i+1])
		val[i] = dat[3][:,i].mean()
	ax[0].plot(rad,val,'b',label='t=10kyr')

	dat = load_2D_dust_data('/data/fung/',304,754,'c50_p21_s5000_r100_i40_o500_ad_rot_fargo_bc10_dr_PEMEU', 1000)
	rad = np.zeros(304)
	val = np.zeros(304)
	for i in range(0,304):
		rad[i] = 0.5*(dat[1][i]+dat[1][i+1])
		val[i] = dat[3][:,i].mean()
	ax[0].plot(rad,val,'r',label='t=100kyr')

	ax[0].legend(fontsize=16)
	ax[0].set_yscale('log')
	ax[0].set_xlabel('radius [au]', fontsize=16)
	ax[0].set_ylabel('surface density', fontsize=16)
	ax[0].tick_params(axis='x', labelsize=16)
	ax[0].tick_params(axis='y', labelsize=16)

	time,ring = gas_ring_2D()
	gas_ring_plot(f1,ax[1],time,ring)
	f1.tight_layout()
	return f1

def gas_dust_3_panels():
	f1,ax=plt.subplots(2, 3, figsize=(14,8), sharex=True, sharey=True)

	dat = load_2D_dust_data('/data/fung/',304,754,'c50_p21_s5000_r100_i40_o500_b10_rot_fargo_PEMEU', 100)
	m1 = make_cart_plot_2D(f1,ax[0,0],dat[0],dat[1],dat[2],dat[3],-1.5,0.5)
	m2 = make_cart_plot_2D(f1,ax[1,0],dat[0],dat[1],dat[2],dat[7]/100.0,-4.0,0.5)

	dat = load_2D_dust_data('/data/fung/',304,754,'c50_p21_s5000_r100_i40_o500_ad_rot_fargo_bc10_PEMEU', 100)
	m1 = make_cart_plot_2D(f1,ax[0,1],dat[0],dat[1],dat[2],dat[3],-1.5,0.5)
	m2 = make_cart_plot_2D(f1,ax[1,1],dat[0],dat[1],dat[2],dat[7]/100.0,-4.0,0.5)

	dat = load_2D_dust_data('/data/fung/',304,754,'c50_p21_s5000_r100_i40_o500_ad_rot_fargo_bc10_dr_PEMEU', 100)
	m1 = make_cart_plot_2D(f1,ax[0,2],dat[0],dat[1],dat[2],dat[3],-1.5,0.5)
	m2 = make_cart_plot_2D(f1,ax[1,2],dat[0],dat[1],dat[2],dat[7]/100.0,-4.0,0.5)

	ax[0,0].set_title('isothermal',fontsize=16)
	ax[0,1].set_title('cooling',fontsize=16)
	ax[0,2].set_title('cooling + dust feedback',fontsize=16)

	ax[0,0].set_xlabel('')
	ax[0,1].set_xlabel('')
	ax[0,2].set_xlabel('')

	ax[0,1].set_ylabel('')
	ax[0,2].set_ylabel('')
	ax[1,1].set_ylabel('')
	ax[1,2].set_ylabel('')

	f1.text(0.0,0.72,'Gas', fontsize=16, rotation=90)
	f1.text(0.0,0.26,'Dust', fontsize=16, rotation=90)

	f1.tight_layout()

	cbar1 = f1.colorbar(m1, ax=ax[0,:].ravel().tolist(), shrink=0.6, location='right', pad=0.02)
	cbar1.ax.tick_params(labelsize=14, width=1)

	cbar2 = f1.colorbar(m2, ax=ax[1,:].ravel().tolist(), shrink=0.6, location='right', pad=0.02)
	cbar2.ax.tick_params(labelsize=14, width=1)

	f1.savefig('../images/gas_dust_compare.png', dpi=300, bbox_inches='tight')
	
	return f1

def res_plot():
	f1,ax=plt.subplots(1, 1, figsize=(6,4.5))

	dat = load_2D_dust_data('/data/fung/',304,754,'c50_p42_s5000_r100_i40_o500_ad_rot_fargo_bc10_dr_PEMEU', 10)
	rad = np.zeros(304)
	val = np.zeros(304)
	for i in range(0,304):
		rad[i] = 0.5*(dat[1][i]+dat[1][i+1])
		val[i] = dat[3][:,i].mean()
	ax.plot(rad,val,'k',label='6 cells/h')

	dat = load_2D_dust_data('/data/fung/',608,1508,'c50_p42_s5000_r100_i40_o500_ad_rot_fargo_bc10_dr_PEMEU', 10)
	rad = np.zeros(608)
	val = np.zeros(608)
	for i in range(0,608):
		rad[i] = 0.5*(dat[1][i]+dat[1][i+1])
		val[i] = dat[3][:,i].mean()
	ax.plot(rad,val,'--b',label='12 cells/h')

	dat = load_2D_dust_data('/data/fung/',1216,3016,'c50_p42_s5000_r100_i40_o500_ad_rot_fargo_bc10_dr_PEMEU', 10)
	rad = np.zeros(1216)
	val = np.zeros(1216)
	for i in range(0,1216):
		rad[i] = 0.5*(dat[1][i]+dat[1][i+1])
		val[i] = dat[3][:,i].mean()
	ax.plot(rad,val,':r',label='24 cells/h')

	ax.legend(fontsize=16)
	ax.set_yscale('log')
	ax.set_ylim(0.15,3.0)
	ax.set_yticks([0.25, 0.5, 1.0, 2.0])
	ax.set_xlabel('radius [au]', fontsize=16)
	ax.set_ylabel('surface density', fontsize=16)
	ax.tick_params(axis='x', labelsize=16)
	ax.tick_params(axis='y', labelsize=16)
	ax.get_yaxis().set_major_formatter(ticker.ScalarFormatter())

	f1.tight_layout()
	return f1
