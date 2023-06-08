import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd
from matplotlib.ticker import ScalarFormatter, NullFormatter  

plt.rcParams['text.usetex'] = True

label = 'h50_PEM3'
xmax = 1024
frame = 1001

def get_growth_rates(m,xmax,label):
	cos_t = np.zeros([m.size,xmax,xmax])
	sin_t = np.zeros([m.size,xmax,xmax])
	dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,xmax,label,0)
	xc = rd.cell_center(dat[1])*np.pi*2.0

	a_min = 1e-10
	a_max = 1e-9

	for n in range(m.size):
		for i in range(xmax):
			sin_t[n,i,:] = np.sin(xc*m[n])
			cos_t[n,i,:] = np.cos(xc*m[n])	


	start = np.zeros(m.size)
	end = np.zeros(m.size)

	pre_amp = np.zeros(m.size)
	amp = np.zeros(m.size)

	amp_start = np.zeros(m.size)
	amp_end = np.zeros(m.size)

	y = np.zeros([m.size,frame])
	x = np.zeros(frame)

	t = 0.0

	dat0 = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,xmax,label,0)

	for i in range(frame):
		pre_t = t
		dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,xmax,label,i)
		vel = np.sqrt(dat[6]*dat[6] + (dat[5]-dat0[5])*(dat[5]-dat0[5]))

		t = dat[0]
		print(i, t)
		x[i] = t

		for n in range(m.size):
			pre_amp[n] = amp[n]

			tmp1 = dat[6]*sin_t[n,:,:]
			tmp2 = dat[6]*cos_t[n,:,:]

			tot_sin = np.sum(tmp1,axis=1)/xmax
			tot_cos = np.sum(tmp2,axis=1)/xmax
			tot = np.sqrt(tot_sin*tot_sin+tot_cos*tot_cos)

			amp[n] = np.sum(tot)/xmax
			y[n,i] = amp[n]

			if (np.max(vel)>1.0e-6 and start[n]==0.0):
				#start[n] = pre_t + (amp[n]-a_min)*(t-pre_t)/(amp[n]-pre_amp[n])
				amp_start[n] = amp[n]
				start[n] = t
	
			if (np.max(vel)>1.0e-5 and end[n]==0.0):
				#end[n] = pre_t + (amp[n]-a_max)*(t-pre_t)/(amp[n]-pre_amp[n])
				amp_end[n] = amp[n]
				end[n] = t

		for n in range(m.size):
			print(y[n,i])

		print(np.min(vel),np.max(vel))

		if (np.max(vel)>1.0e-2):
			break

	for n in range(m.size):
		print(start[n],end[n],amp_start[n],amp_end[n],(np.log(amp_end[n]/amp_start[n]))/(end[n]-start[n]))

	for n in range(8):
		plt.plot(x,y[n,:],label='m='+str(m[n]))
	plt.yscale('log')
	plt.legend()
	plt.show()

#m=np.array([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0])#,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0])
#get_growth_rates(m,xmax,label)

#dat0 = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,xmax,label,0)
dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,xmax,label,65)
plt.pcolormesh(dat[1],dat[2],dat[3])
plt.show()

#for i in range(0,frame+1):
#	dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax[0],xmax[0],label_pp4,i)

#	fig, axs = plt.subplots(1, 1, figsize=(8,8))
#	axs.pcolormesh(dat[1],dat[2],dat[3],cmap='gnuplot',vmin=0.9,vmax=2.1)

#	fname = '/home/fung/PEnGUIn/images/PP4_'+str(xmax[0])+'_'+rd.frame_num(i)+'.png'
#	fig.savefig(fname, bbox_inches='tight', dpi=600)
#	print(fname)
#	plt.close()
