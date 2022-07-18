import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd

xmax = 480
ymax = 480
label = 'h50'

pem_v_ener = np.zeros(501)
plm_v_ener = np.zeros(501)
ppm_v_ener = np.zeros(501)
sin_template = np.zeros((xmax,ymax));
cos_template = np.zeros((xmax,ymax));
u_template = np.zeros((xmax,ymax));

dat = rd.load_2D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,label,0)
xc = rd.cell_center(dat[1])
yc = rd.cell_center(dat[2])

for i in range(0,ymax):
	sin_template[i,:] = np.sin(16.0*np.pi*xc)
	cos_template[i,:] = np.cos(16.0*np.pi*xc)

for i in range(0,ymax):
	if (yc[i]>=0.25 and yc[i]<=0.75): u_template[i,:] = 0.5
	else: u_template[i,:] = -0.5

#plt.figure(1)
#plt.pcolormesh(dat[1],dat[2],dat[5]-u_template);

#plt.figure(2)
#plt.pcolormesh(dat[1],dat[2],dat[6]);
#plt.show()


for i in range(0,501):
	frame = i
	plm = rd.load_2D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,label+'_PLM',frame)
	pem = rd.load_2D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,label,frame)
	ppm = rd.load_2D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,label+'_PPM',frame)

	stmp = pem[3]*pem[6]*sin_template;
	ctmp = pem[3]*pem[6]*cos_template;
	pem_v_ener[i] = np.sum(stmp*stmp + ctmp*ctmp)/xmax/ymax

	stmp = plm[3]*plm[6]*sin_template;
	ctmp = plm[3]*plm[6]*cos_template;
	plm_v_ener[i] = np.sum(stmp*stmp + ctmp*ctmp)/xmax/ymax

	stmp = ppm[3]*ppm[6]*sin_template;
	ctmp = ppm[3]*ppm[6]*cos_template;
	ppm_v_ener[i] = np.sum(stmp*stmp + ctmp*ctmp)/xmax/ymax

	print(i, plm_v_ener[i],pem_v_ener[i],ppm_v_ener[i])

#	fig, axs = plt.subplots(1, 3, figsize=(36,14), sharey=True)

#	pcm = axs[0].pcolormesh(plm[1],plm[2],plm[3],cmap='gnuplot',vmin=0.9,vmax=2.1)
#	axs[0].set_aspect('equal')
#	axs[0].tick_params(axis='x', labelsize=16)
#	axs[0].tick_params(axis='y', labelsize=16)

#	pcm = axs[1].pcolormesh(pem[1],pem[2],pem[3],cmap='gnuplot',vmin=0.9,vmax=2.1)
#	axs[1].set_aspect('equal')
#	axs[1].tick_params(axis='x', labelsize=16)

#	pcm = axs[2].pcolormesh(ppm[1],ppm[2],ppm[3],cmap='gnuplot',vmin=0.9,vmax=2.1)
#	axs[2].set_aspect('equal')
#	axs[2].tick_params(axis='x', labelsize=16)
	
#	fig.subplots_adjust(wspace=0,hspace=0)

#	cbar = fig.colorbar(pcm, ax=axs[:], shrink=0.4, pad=0.05, location='bottom')
#	cbar.ax.tick_params(labelsize=16)

#	fig.savefig('../images/frame_KH_mres_'+rd.frame_num(i)+'.png', bbox_inches='tight')

#	plt.close()
#	print(i)

plt.plot(np.sqrt(pem_v_ener),'r')
plt.plot(np.sqrt(plm_v_ener),'b')
plt.plot(np.sqrt(ppm_v_ener),'g')
plt.yscale('log')

plt.show()
