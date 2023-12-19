import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
from matplotlib import cm
import numpy as np
import math
import os.path as pat

twopi = 2.0*np.pi

def check_file_exist_2D(path, imax, jmax, label, num):
	fname = path+'binary_'
	fname = fname + str(imax)+'x'+str(jmax)+'_'
	fname = fname + label + '_'
	fname = fname + frame_num(num)
	return pat.isfile(fname)

def check_file_exist_3D(path, imax, jmax, kmax, label, num):
	fname = path+'binary_'
	fname = fname + str(imax)+'x'+str(jmax)+'x'+str(kmax)+'_'
	fname = fname + label + '_'
	fname = fname + frame_num(num)
	return pat.isfile(fname)

def error(a,n):
	if (n==-1):
		return np.max(np.abs(a))
	elif (n==1):
		return np.sum(np.abs(a))/a.size
	else:
		return np.power(np.mean(np.power(a,n)),1.0/n)

def cell_center(xa):
	N = len(xa)-1
	xc = np.arange(N,dtype=np.float64)

	for i in range(N):
		xc[i] = (xa[i+1]+xa[i])/2.0
	return xc

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
	x = dat[start:end]

	start = end
	end = start+imax
	ro = dat[start:end]

	start = end
	end = start+imax
	pr = dat[start:end]

	start = end
	end = start+imax
	vx = dat[start:end]

	return (t,x,ro,pr,vx)

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

def load_1D_dust_data(path, imax, label, num):
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

	start = end
	end = start+imax
	rd = dat[start:end]

	start = end
	end = start+imax
	xd = dat[start:end]

	return (t,xc,ro,pr,vx,rd,xd)

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
	rd = np.reshape(dat[start:end], (-1,imax))

	start = end
	end = start+imax*jmax
	xd = np.reshape(dat[start:end], (-1,imax))

	start = end
	end = start+imax*jmax
	yd = np.reshape(dat[start:end], (-1,imax))

	return (time,x,y,ro,pr,vx,vy,rd,xd,yd)
	
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
	rd = np.reshape(dat[start:end], (-1,jmax,imax))

	start = end
	end = start+imax*jmax*kmax
	xd = np.reshape(dat[start:end], (-1,jmax,imax))

	start = end
	end = start+imax*jmax*kmax
	yd = np.reshape(dat[start:end], (-1,jmax,imax))
	
	start = end
	end = start+imax*jmax*kmax
	zd = np.reshape(dat[start:end], (-1,jmax,imax))

	return (time,x,y,z,ro,pr,vx,vy,vz,rd,xd,yd,zd)

def make_cart_plot_2D(fig,ax,time,x,y,val,vmin,vmax):
	pres = 512
	imax = len(x)-1
	jmax = len(y)-1
	
	rad = np.zeros(imax)
	azi = np.zeros(jmax)
    
	origin = [-1.0,0.0]
	box_size = 0.1
	cx = np.zeros(pres+1)
	cy = np.zeros(pres+1)
	dd = 0.5*box_size/pres
	aa = np.reshape(np.zeros(pres*pres), (-1,pres))
	
	for i in range(0,pres+1):
		cx[i] = (box_size*i)/pres - 0.5*box_size + origin[0];
		cy[i] = (box_size*i)/pres - 0.5*box_size + origin[1];

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
	ax.set_xticks(np.arange(-box_size + origin[0], box_size + origin[0], step=box_size/4))
	ax.set_yticks(np.arange(-box_size + origin[1], box_size + origin[1], step=box_size/4))
	ax.set_xlim(-0.5*box_size + origin[0],0.5*box_size + origin[0])
	ax.set_ylim(-0.5*box_size + origin[1],0.5*box_size + origin[1])
	ax.tick_params(axis='both', labelsize=12, width=2)

	return mesh1

def get_val(dat):
	return dat[4]/dat[3]    

def make_figure_2D(p0,p1,p2,p3):
	val0 = get_val(p0)
	val1 = get_val(p1)
	val2 = get_val(p2)
	val3 = get_val(p3)
    
	fig, axs = plt.subplots(2, 2, figsize=(11,11))
	vmin = np.min([np.min(val0),np.min(val1),np.min(val2),np.min(val3)])
	vmax = np.max([np.max(val0),np.max(val1),np.max(val2),np.max(val3)])
	vmin = np.log10(vmin)
	vmax = np.log10(vmax)

	m1=make_cart_plot_2D(fig,axs[0,0],p0[0],p0[1],p0[2],val0,vmin,vmax)
	m2=make_cart_plot_2D(fig,axs[0,1],p1[0],p1[1],p1[2],val1,vmin,vmax)
	m3=make_cart_plot_2D(fig,axs[1,0],p2[0],p2[1],p2[2],val2,vmin,vmax)
	m4=make_cart_plot_2D(fig,axs[1,1],p3[0],p3[1],p3[2],val3,vmin,vmax)

	axs[0,0].set_title('1.')
	axs[0,1].set_title('2.')
	axs[1,0].set_title('3.')
	axs[1,1].set_title('4.')

	#plt.figtext(0.41,0.88,'t = '+str(round(p0[0]/(2.0*np.pi),2)/1000.0)+' kyr', fontsize=18)
	cbar1 = fig.colorbar(m1, ax=axs[0,0], shrink=0.6, pad=0.02)
	cbar1.ax.tick_params(labelsize=12, width=1)

	cbar2 = fig.colorbar(m2, ax=axs[0,1], shrink=0.6, pad=0.02)
	cbar2.ax.tick_params(labelsize=12, width=1)
	
	cbar2 = fig.colorbar(m3, ax=axs[1,0], shrink=0.6, pad=0.02)
	cbar2.ax.tick_params(labelsize=12, width=1)
	
	cbar2 = fig.colorbar(m4, ax=axs[1,1], shrink=0.6, pad=0.02)
	cbar2.ax.tick_params(labelsize=12, width=1)
	return fig,axs
