import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd
<<<<<<< HEAD

label_pem = 'h50'
label_plm = 'h50_PLM'
label_ppm = 'h50_PPM'
xmax = 128

frame = 10

#std = rd.load_1D_data('/mnt/penguin/fung/p2/',3200,label_pem,frame)
#hpe = rd.load_1D_data('/mnt/penguin/fung/p2/',400,label_pem,frame)
#hpl = rd.load_1D_data('/mnt/penguin/fung/p2/',400,label_plm,frame)
#hpp = rd.load_1D_data('/mnt/penguin/fung/p2/',400,label_ppm,frame)
pem = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_pem,frame)
plm = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_plm,frame)
ppm = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_ppm,frame)

#pem_r = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_pem+'_rev',frame)
#plm_r = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_plm+'_rev',frame)
#ppm_r = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_ppm+'_rev',frame)
#a = np.flip(pem_r[2])
#b = np.flip(plm_r[2])
#c = np.flip(ppm_r[2])

xc = pem[1]
sol = 1.0 + 1.e-6*np.sin(2.0*np.pi*xc)

plt.figure(1)

#plt.plot(std[1],std[2],'grey')


plt.plot(xc,ppm[2]-sol,'b:', mfc='none')
#plt.plot(hpp[1],hpp[2],'bx-')
#plt.plot(xc,c,'bx')
#print(np.max(np.abs(ppm[2]-c)))

plt.plot(xc,plm[2]-sol,'g:', mfc='none')
#plt.plot(hpl[1],hpl[2],'gx-')
#plt.plot(xc,b,'gx')
#print(np.max(np.abs(plm[2]-b)))

plt.plot(xc,pem[2]-sol,'r:', mfc='none')
#plt.plot(hpe[1],hpe[2],'rx-')
#plt.plot(xc,a,'rx')
#print(np.max(np.abs(pem[2]-a)))

pem = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_pem,frame)
plm = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_plm,frame)
ppm = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_ppm,frame)

xc = pem[1]
plt.figure(2)

plt.plot(xc,ppm[3],'bs-', mfc='none')

plt.plot(xc,plm[3],'gs-', mfc='none')

plt.plot(xc,pem[3],'rs-', mfc='none')

pem = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_pem,frame)
plm = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_plm,frame)
ppm = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_ppm,frame)

xc = pem[1]

plt.figure(3)

plt.plot(xc,ppm[4],'bs-', mfc='none')

plt.plot(xc,plm[4],'gs-', mfc='none')

plt.plot(xc,pem[4],'rs-', mfc='none')


plt.show()
=======
from matplotlib.ticker import ScalarFormatter, NullFormatter  

plt.rcParams['text.usetex'] = True

CFL = np.array(['_CFL5','_CFL20','_CFL80'])
CFL_val = np.array([0.05,0.2,0.8])

label_moc = 'h50_MOC'
label_van = 'h50_VAN'
label_pe3 = 'h50_PEM3'
label_pp3 = 'h50_PPM3'
label_pe4 = 'h50_PEM4'
label_pp4 = 'h50_PPM4'

xmax = np.array([32,64,128,256])
frame = 1

a0 = np.zeros([4,4])
a1 = np.zeros([4,4])
a2 = np.zeros([4,4])
a3 = np.zeros([4,4])
a4 = np.zeros([4,4])
a5 = np.zeros([4,4])

order = 2

for j in range(xmax.size):
	van = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax[j],label_van+CFL[0],0)
	xc = van[1]
	sol = 1.0 + 1.e-6*np.sin(2.0*np.pi*xc)

	for i in range(CFL_val.size):
		van = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax[j],label_van+CFL[i],frame)
		moc = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax[j],label_moc+CFL[i],frame)
		pe3 = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax[j],label_pe3+CFL[i],frame)
		pp3 = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax[j],label_pp3+CFL[i],frame)
		pe4 = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax[j],label_pe4+CFL[i],frame)
		pp4 = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax[j],label_pp4+CFL[i],frame)
	
		a0[i,j] = rd.error(van[2]-sol,order)
		a1[i,j] = rd.error(moc[2]-sol,order)
		a2[i,j] = rd.error(pp3[2]-sol,order)
		a3[i,j] = rd.error(pe3[2]-sol,order)
		a4[i,j] = rd.error(pp4[2]-sol,order)
		a5[i,j] = rd.error(pe4[2]-sol,order)

style = np.array(['b+--','g+--','r^--','k^-','rs:','ks:'])
label = np.array(['VAN2','MOC2','PPM3','PEM3','PPM4','PEM4'])

fig, axs = plt.subplots(1, 3, figsize=(16*0.9,9*0.45), sharey='row')

for i in range(3):
	axs[i].plot(xmax,a0[i,:],style[0], mfc='none',label=label[0])
	axs[i].plot(xmax,a1[i,:],style[1], mfc='none',label=label[1])
	axs[i].plot(xmax,a2[i,:],style[2], mfc='none',label=label[2])
	axs[i].plot(xmax,a3[i,:],style[3], mfc='none',label=label[3])
	axs[i].plot(xmax,a4[i,:],style[4], mfc='none',label=label[4])
	axs[i].plot(xmax,a5[i,:],style[5], mfc='none',label=label[5])
	axs[i].plot(xmax[0:3],5e-10*xmax[0]*xmax[0]/xmax[0:3]/xmax[0:3],color='grey',linestyle='--')
	axs[i].legend()
	axs[i].set_xlabel(r'\# of Cells')
	
	tt = r'\textbf{CFL Number = '+str(CFL_val[i])+'}'
	axs[i].text(0.05,0.1,tt, ha='left', va='center', transform=axs[i].transAxes)

	tt = r'$2^{\rm nd}$ order'
	axs[i].text(0.21,0.28,tt, ha='left', va='center', transform=axs[i].transAxes, rotation=-25)
		

#	axs[1,i].plot(CFL_val,a0[:,i],style[0], mfc='none',label=label[0])
#	axs[1,i].plot(CFL_val,a1[:,i],style[1], mfc='none',label=label[1])
#	axs[1,i].plot(CFL_val,a2[:,i],style[2], mfc='none',label=label[2])
#	axs[1,i].plot(CFL_val,a3[:,i],style[3], mfc='none',label=label[3])
#	axs[1,i].plot(CFL_val,a4[:,i],style[4], mfc='none',label=label[4])
#	axs[1,i].plot(CFL_val,a5[:,i],style[5], mfc='none',label=label[5])
#	axs[1,i].legend()
#	axs[1,i].set_xlabel('CFL Number')

#	tt = r'\textbf{\# of Cells = '+str(xmax[i])+'}'
#	axs[1,i].text(0.05,0.9,tt, ha='left', va='center', transform=axs[1,i].transAxes)

if (order==-1):
	axs[0].set_ylabel(r'\textbf{L{$\infty$} Error}')
#	axs[1,0].set_ylabel(r'\textbf{L{$\infty$} Error}')
else:
	axs[0].set_ylabel(r'\textbf{L'+str(order)+' Error}')
#	axs[1,0].set_ylabel(r'\textbf{L'+str(order)+' Error}')

for i in range(3):
	axs[i].set_xscale('log')
	for axis in [axs[i].xaxis, axs[i].yaxis]:
		axis.set_major_formatter(ScalarFormatter())
		axis.set_minor_formatter(NullFormatter())
	axs[i].set_xticks(xmax)
	axs[i].set_yscale('log')

#	axs[1,i].set_xscale('log')
#	for axis in [axs[1,i].xaxis, axs[1,i].yaxis]:
#		axis.set_major_formatter(ScalarFormatter())
#		axis.set_minor_formatter(NullFormatter())
#	axs[1,i].set_xticks(CFL_val)
#	axs[1,i].set_yscale('log')

fig.savefig('L2_linear.png', bbox_inches='tight', dpi=600)
>>>>>>> methods
