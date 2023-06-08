import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd
from matplotlib.ticker import ScalarFormatter, NullFormatter  

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 16})

label = ['h35_3p10E_OA_MOC_2dev','h35_4p10E_OA_MOC_2dev','h35_5p10E_OA_MOC_2dev','h35_6p10E_OA_MOC_2dev','h35_7p10E_OA_MOC_2dev','h35_1p100E_OA_MOC_2dev']
legend= ['3p','4p','5p','6p','7p','1p']

ndat = 6
dat = [None]*ndat
xmax = 960
ymax = 1632
start = 1000

#label = ['h100_1p95E_a-40_OA_MOC','h100_1p95E_OA_MOC']
#legend= [r'$\alpha$=0.0001','inviscid']

#ndat = 2
#dat = [None]*ndat
#xmax = 456
ymax = 1008
start = 4000


plt.figure(10)
rad = np.logspace(np.log10(0.3),np.log10(10.0),xmax)
#plt.plot(rad,np.power(rad,-1.5),"k--")
plt.ylabel(r'$\Delta\Sigma_{\rm gas}$')
plt.xlabel('radius [au]')

def circle(ax,origin,radius,style):
	angles = np.linspace(0.0,2.0*np.pi,150)
	x = origin[0]+radius*np.cos(angles)
	y = origin[1]+radius*np.sin(angles)
	ax.plot(x,y,style,linewidth=1.0)
	return

def cartesian(dat,smax,res):
	xa = np.linspace(-smax,smax,res)
	ya = np.linspace(-smax,smax,res)
	xc = rd.cell_center(xa)
	yc = rd.cell_center(ya)

	rad = np.zeros([xc.size,yc.size])
	phi = np.zeros([xc.size,yc.size])
	rho = np.zeros([xc.size,yc.size])

	rc = rd.cell_center(dat[1])
	pc = rd.cell_center(dat[2])

	for i in range(xc.size):
		rad[i,:] = np.sqrt(xc[i]*xc[i]+yc*yc)
		phi[i,:] = np.arctan2(yc,xc[i])
		for j in range(yc.size):
			if rad[i,j]>0.7:
				rho[i,j] = rd.interpolate_2D_periodicY([rad[i,j],phi[i,j]],rc,pc,dat[3],0.0,2.0*np.pi)
			else:
				rho[i,j] = None			
			#if rho[i,j] is not None:
			#	rho[i,j]*=np.power(rad[i,j],1.5)
	return xa,ya,rho
	

for i in range(0,ndat):
	plt.figure(i)

	dat[i] = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,ymax,label[i],start)
	xa,ya,rho = cartesian(dat[i],3.0,512)
	#plt.pcolormesh(dat[i][1],dat[i][2],np.log10(dat[i][3]))
	plt.pcolormesh(36.0*xa,36.0*ya,rho,vmin=0.1,vmax=1.6,cmap='pink')

	ax = plt.gca()
	circle(ax,[0,0],36.0,'w--')
	circle(ax,[0,0],60.0,'g--')
	ax.set_aspect('equal')

	plt.colorbar()

	plt.ylabel('y [au]')
	plt.xlabel('x [au]')
	plt.title(legend[i])

	plt.figure(10)
	rad = rd.cell_center(dat[i][1])
	plt.plot(36.0*rad,np.mean(dat[i][3],axis=0),label=legend[i])

	#plt.plot(rad,np.mean(dat[3],axis=0)/np.power(rad,-1.5)-1,label=legend[i])

#plt.figure(11)
#plt.plot(np.mean(dat[5][3]-dat[6][3],axis=0))

plt.figure(10)
plt.ylim([0.0,1.0])
plt.xlim([0.0,300.0])
plt.legend()
plt.show()

#for i in range(1,1844,4):
#	dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,ymax,label,i)
#	print(i,dat[0],np.min(np.mean(dat[3],axis=0)))
	
