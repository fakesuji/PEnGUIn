import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd
import matplotlib.ticker

plt.rcParams['text.usetex'] = True

CFL = np.array(['_CFL10','_CFL20','_CFL40','_CFL80'])
CFL_val = np.array([0.1,0.2,0.4,0.8])

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

for j in range(4):
	van = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax[j],label_van+CFL[0],0)
	xc = van[1]
	sol = 1.0 + 1.e-6*np.sin(2.0*np.pi*xc)

	for i in range(4):
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

fig, axs = plt.subplots(2, 4, figsize=(24*0.9,9*0.9), sharey='row')

for i in range(4):
	axs[0,i].plot(xmax,a0[i,:],style[0], mfc='none',label=label[0])
	axs[0,i].plot(xmax,a1[i,:],style[1], mfc='none',label=label[1])
	axs[0,i].plot(xmax,a2[i,:],style[2], mfc='none',label=label[2])
	axs[0,i].plot(xmax,a3[i,:],style[3], mfc='none',label=label[3])
	axs[0,i].plot(xmax,a4[i,:],style[4], mfc='none',label=label[4])
	axs[0,i].plot(xmax,a5[i,:],style[5], mfc='none',label=label[5])
	axs[0,i].legend()
	axs[0,i].set_xlabel(r'\# of Cells')
	
	tt = r'\textbf{CFL Number = '+str(CFL_val[i])+'}'
	axs[0,i].text(0.05,0.1,tt, ha='left', va='center', transform=axs[0,i].transAxes)
		

	axs[1,i].plot(CFL_val,a0[:,i],style[0], mfc='none',label=label[0])
	axs[1,i].plot(CFL_val,a1[:,i],style[1], mfc='none',label=label[1])
	axs[1,i].plot(CFL_val,a2[:,i],style[2], mfc='none',label=label[2])
	axs[1,i].plot(CFL_val,a3[:,i],style[3], mfc='none',label=label[3])
	axs[1,i].plot(CFL_val,a4[:,i],style[4], mfc='none',label=label[4])
	axs[1,i].plot(CFL_val,a5[:,i],style[5], mfc='none',label=label[5])
	axs[1,i].legend()
	axs[1,i].set_xlabel('CFL Number')

	tt = r'\textbf{\# of Cells = '+str(xmax[i])+'}'
	axs[1,i].text(0.05,0.9,tt, ha='left', va='center', transform=axs[1,i].transAxes)

if (order==-1):
	axs[0,0].set_ylabel(r'\textbf{L{$\infty$} Error}')
	axs[1,0].set_ylabel(r'\textbf{L{$\infty$} Error}')
else:
	axs[0,0].set_ylabel(r'\textbf{L'+str(order)+' Error}')
	axs[1,0].set_ylabel(r'\textbf{L'+str(order)+' Error}')

for i in range(4):
	axs[0,i].set_yscale('log')
	#axs[0,i].set_xscale('log')
	axs[0,i].set_xticks(xmax)
	axs[0,i].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

	axs[1,i].set_yscale('log')
	#axs[1,i].set_xscale('log')
	axs[1,i].set_xticks(CFL_val)
	axs[1,i].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

fig.savefig('tmp.png', bbox_inches='tight')
