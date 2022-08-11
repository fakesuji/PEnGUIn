import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd
from matplotlib.ticker import ScalarFormatter, NullFormatter  

plt.rcParams['text.usetex'] = True

label = 'h50_PEM3'
xmax = 1024
frame = 700

t = np.zeros(frame)
amp = np.zeros(frame)

cos_t = np.zeros([xmax,xmax])
sin_t = np.zeros([xmax,xmax])
dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,xmax,label,0)
xc = rd.cell_center(dat[1])*np.pi*2.0

for i in range(xmax):
	sin_t[i,:] = np.sin(xc)
	cos_t[i,:] = np.cos(xc)

for i in range(frame):
	dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,xmax,label,i+1)
	tmp1 = 2.0*np.sum(dat[3]*sin_t)/xmax/xmax
	tmp2 = 2.0*np.sum(dat[3]*cos_t)/xmax/xmax
	t[i] = dat[0]
	amp[i] = np.sqrt(tmp1*tmp1+tmp2*tmp2)
	print(t[i],amp[i])

plt.plot(t, amp, label='m=1')

for i in range(xmax):
	sin_t[i,:] = np.sin(xc*2.0)
	cos_t[i,:] = np.cos(xc*2.0)

for i in range(frame):
	dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,xmax,label,i+1)
	tmp1 = 2.0*np.sum(dat[3]*sin_t)/xmax/xmax
	tmp2 = 2.0*np.sum(dat[3]*cos_t)/xmax/xmax
	t[i] = dat[0]
	amp[i] = np.sqrt(tmp1*tmp1+tmp2*tmp2)
	print(t[i],amp[i])

plt.plot(t, amp, label='m=2')

for i in range(xmax):
	sin_t[i,:] = np.sin(xc*3.0)
	cos_t[i,:] = np.cos(xc*3.0)

for i in range(frame):
	dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,xmax,label,i+1)
	tmp1 = 2.0*np.sum(dat[3]*sin_t)/xmax/xmax
	tmp2 = 2.0*np.sum(dat[3]*cos_t)/xmax/xmax
	t[i] = dat[0]
	amp[i] = np.sqrt(tmp1*tmp1+tmp2*tmp2)
	print(t[i],amp[i])

plt.plot(t, amp, label='m=3')

for i in range(xmax):
	sin_t[i,:] = np.sin(xc*4.0)
	cos_t[i,:] = np.cos(xc*4.0)

for i in range(frame):
	dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,xmax,label,i+1)
	tmp1 = 2.0*np.sum(dat[3]*sin_t)/xmax/xmax
	tmp2 = 2.0*np.sum(dat[3]*cos_t)/xmax/xmax
	t[i] = dat[0]
	amp[i] = np.sqrt(tmp1*tmp1+tmp2*tmp2)
	print(t[i],amp[i])

plt.plot(t, amp, label='m=4')

for i in range(xmax):
	sin_t[i,:] = np.sin(xc*5.0)
	cos_t[i,:] = np.cos(xc*5.0)

for i in range(frame):
	dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,xmax,label,i+1)
	tmp1 = 2.0*np.sum(dat[3]*sin_t)/xmax/xmax
	tmp2 = 2.0*np.sum(dat[3]*cos_t)/xmax/xmax
	t[i] = dat[0]
	amp[i] = np.sqrt(tmp1*tmp1+tmp2*tmp2)
	print(t[i],amp[i])

plt.plot(t, amp, label='m=5')

for i in range(xmax):
	sin_t[i,:] = np.sin(xc*6.0)
	cos_t[i,:] = np.cos(xc*6.0)

for i in range(frame):
	dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,xmax,label,i+1)
	tmp1 = 2.0*np.sum(dat[3]*sin_t)/xmax/xmax
	tmp2 = 2.0*np.sum(dat[3]*cos_t)/xmax/xmax
	t[i] = dat[0]
	amp[i] = np.sqrt(tmp1*tmp1+tmp2*tmp2)
	print(t[i],amp[i])

plt.plot(t, amp, label='m=6')

for i in range(xmax):
	sin_t[i,:] = np.sin(xc*7.0)
	cos_t[i,:] = np.cos(xc*7.0)

for i in range(frame):
	dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,xmax,label,i+1)
	tmp1 = 2.0*np.sum(dat[3]*sin_t)/xmax/xmax
	tmp2 = 2.0*np.sum(dat[3]*cos_t)/xmax/xmax
	t[i] = dat[0]
	amp[i] = np.sqrt(tmp1*tmp1+tmp2*tmp2)
	print(t[i],amp[i])

plt.plot(t, amp, label='m=7')

for i in range(xmax):
	sin_t[i,:] = np.sin(xc*8.0)
	cos_t[i,:] = np.cos(xc*8.0)

for i in range(frame):
	dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,xmax,label,i+1)
	tmp1 = 2.0*np.sum(dat[3]*sin_t)/xmax/xmax
	tmp2 = 2.0*np.sum(dat[3]*cos_t)/xmax/xmax
	t[i] = dat[0]
	amp[i] = np.sqrt(tmp1*tmp1+tmp2*tmp2)
	print(t[i],amp[i])

plt.plot(t, amp, label='m=8')

plt.yscale('log')
plt.legend()
plt.show()

#for i in range(0,frame+1):
#	dat = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax[0],xmax[0],label_pp4,i)

#	fig, axs = plt.subplots(1, 1, figsize=(8,8))
#	axs.pcolormesh(dat[1],dat[2],dat[3],cmap='gnuplot',vmin=0.9,vmax=2.1)

#	fname = '/home/fung/PEnGUIn/images/PP4_'+str(xmax[0])+'_'+rd.frame_num(i)+'.png'
#	fig.savefig(fname, bbox_inches='tight', dpi=600)
#	print(fname)
#	plt.close()
