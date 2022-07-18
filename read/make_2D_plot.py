import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd

second = 0

label = 'h50_p4J_a-40_OA_St-30'
labeh = 'h50_PLM'
xmax = 576
ymax = 1024

#xmax = 720
#ymax = 1536
#label = 'h25_p12J_a-10_OA_PPM'

#label = 'h50_PPM';#_p3J_a-2_OA'
#labeh = 'h50_PLM'
#xmax = 960
#ymax = 960

frame = 20
inc = 2
p0 = rd.load_2D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,label,frame)
p1 = rd.load_2D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,label,frame//inc)
p2 = rd.load_2D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,label,frame//inc//inc)
p3 = rd.load_2D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,label,frame//inc//inc//inc)
dat0 = rd.load_2D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,label,0)

xc = rd.cell_center(dat0[1])

plt.figure(1)
plt.pcolormesh(dat0[1],dat0[2],np.log10(p0[3]))
plt.colorbar()

plt.figure(2)
plt.plot(xc,(dat0[5]).mean(axis=0),':k')
plt.plot(xc,(p0[5]).mean(axis=0),'r')

plt.figure(3)
plt.pcolormesh(dat0[1],dat0[2],np.log10(p0[7]))
plt.colorbar()

plt.figure(4)
plt.plot(xc,(dat0[8]).mean(axis=0),':k')
plt.plot(xc,(p0[8]).mean(axis=0),'r')

#plt.figure(2)
#plt.pcolormesh(dat0[1],dat0[2],p0[4]/p0[3])
#plt.colorbar()

#plt.figure(3)
#plt.plot(xc,(dat0[3]).mean(axis=0),':k')
#plt.plot(xc,(p3[3]).mean(axis=0),'r')
#plt.plot(xc,(p2[3]).mean(axis=0),'orange')
#plt.plot(xc,(p1[3]).mean(axis=0),'green')
#plt.plot(xc,(p0[3]).mean(axis=0),'blue')
#plt.xscale("log")
#plt.yscale("log")

#plt.figure(4)
#plt.plot(xc,(dat0[5]).mean(axis=0),'k:')
#plt.plot(xc,(p3[5]).mean(axis=0),'r')
#plt.plot(xc,(p2[5]).mean(axis=0),'orange')
#plt.plot(xc,(p1[5]).mean(axis=0),'g')
#plt.plot(xc,(p0[5]).mean(axis=0),'b')

if (second==1):
	#frame = 2400
	s0 = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,ymax,labeh,frame)
	s1 = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,ymax,labeh,frame//inc)
	s2 = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,ymax,labeh,frame//inc//inc)
	s3 = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,ymax,labeh,frame//inc//inc//inc)
	dat0 = rd.load_2D_data('/mnt/penguin/fung/p2/',xmax,ymax,labeh,0)
	xc = rd.cell_center(dat0[1])

	plt.figure(5)
	plt.pcolormesh(dat0[1],dat0[2],s0[3],vmin=vmin, vmax=vmax)
	plt.colorbar()

	plt.figure(6)
	plt.pcolormesh(dat0[1],dat0[2],s0[4]/s0[3])
	plt.colorbar()

	plt.figure(7)
	plt.plot(xc,(dat0[3]).mean(axis=0),':k')
	plt.plot(xc,(s3[3]).mean(axis=0),'r')
	plt.plot(xc,(s2[3]).mean(axis=0),'orange')
	plt.plot(xc,(s1[3]).mean(axis=0),'green')
	plt.plot(xc,(s0[3]).mean(axis=0),'blue')
	plt.xscale("log")
	plt.yscale("log")

	plt.figure(8)
	plt.plot(xc,(dat0[5]).mean(axis=0),'k:')
	plt.plot(xc,(s3[5]).mean(axis=0),'r')
	plt.plot(xc,(s2[5]).mean(axis=0),'orange')
	plt.plot(xc,(s1[5]).mean(axis=0),'g')
	plt.plot(xc,(s0[5]).mean(axis=0),'b')

	comp = s0[3]/p0[3]-1.0

	plt.figure(9)
	plt.plot(xc,(p0[3]).mean(axis=0),'k')
	plt.plot(xc,(s0[3]).mean(axis=0),'blue')
	plt.xscale("log")
	plt.yscale("log")

	#plt.figure(10)
	#plt.plot(comp.mean(axis=0),'r')
	#print(np.min(comp))
	#print(np.max(comp))
else:
	comp = p0[3]/dat0[3]-1.0

	#plt.figure(9)
	#plt.plot(comp.mean(axis=1),'r')

	#plt.figure(10)
	#plt.plot(comp.mean(axis=0),'r')
	#print(np.min(comp))
	#print(np.max(comp))

	#plt.figure(11)
	#plt.pcolormesh(dat0[1],dat0[2],comp)
	#plt.colorbar()

plt.show()
