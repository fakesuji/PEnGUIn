import read as pg
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

#d0 = pg.load_3D_data("/mnt/penguin/fung/p2/",432,1008,32,"h50_p11J_a-2_OA",0)
#d1 = pg.load_3D_data("/mnt/penguin/fung/p2/",432,1008,32,"h50_p11J_a-2_OA",1)
#m1 = pg.load_3D_data("/mnt/penguin/fung/p2/",432,1008,32,"h50_p11J_a-2_OA_multi",1)
#print(np.max(np.abs(d1[4]-m1[4])))
#print(np.max(np.abs(d1[4]-d0[4])))

#plt.figure(1)
#plt.pcolormesh(d1[1],d1[3],d1[4].mean(axis=1)/d0[4].mean(axis=1)-1)
#plt.colorbar()
#plt.figure(2)
#plt.pcolormesh(d1[1],d1[2],d1[4].mean(axis=0)/d0[4].mean(axis=0)-1)
#plt.colorbar()
#plt.figure(3)
#plt.pcolormesh(d1[2],d1[3],d1[4].mean(axis=2)/d0[4].mean(axis=2)-1)
#plt.colorbar()

#plt.figure(4)
#plt.pcolormesh(d1[1],d1[3],m1[4].mean(axis=1)/d0[4].mean(axis=1)-1)
#plt.colorbar()
#plt.figure(5)
#plt.pcolormesh(d1[1],d1[2],m1[4].mean(axis=0)/d0[4].mean(axis=0)-1)
#plt.colorbar()
#plt.figure(6)
#plt.pcolormesh(d1[2],d1[3],m1[4].mean(axis=2)/d0[4].mean(axis=2)-1)
#plt.colorbar()

#plt.figure(7)
#plt.pcolormesh(d1[1],d1[3],m1[4].mean(axis=1)-d1[4].mean(axis=1))
#plt.colorbar()
#plt.figure(8)
#plt.pcolormesh(d1[1],d1[2],m1[4].mean(axis=0)-d1[4].mean(axis=0))
#plt.colorbar()
#plt.figure(9)
#plt.pcolormesh(d1[2],d1[3],m1[4].mean(axis=2)-d1[4].mean(axis=2))
#plt.colorbar()

d0 = pg.load_3D_data("/mnt/penguin/fung/p2/",432,1008,32,"h50_p1J_a-3_OA_multi",0)
d1 = pg.load_3D_data("/mnt/penguin/fung/p2/",432,1008,32,"h50_p1J_a-3_OA",21)
tmp = d1[4]/d0[4]

x0 = 1.0-0.4;
x1 = 1.0+0.4;

y0 = np.pi-0.5*(x1-x0);
y1 = np.pi+0.5*(x1-x0);

plt.figure(1)
plt.pcolormesh(d0[1],d0[2],tmp[31,:,:],norm=matplotlib.colors.LogNorm(),cmap="gnuplot2")
plt.ylim([y0,y1])
plt.xlim([x0,x1])
ax = plt.gca();
ax.set_aspect('equal')
plt.colorbar()

plt.figure(2)
plt.pcolormesh(d0[1],d0[3],d1[4][:,504,:],norm=matplotlib.colors.LogNorm(),cmap="gnuplot2")
plt.colorbar()

plt.figure(3)
plt.pcolormesh(d0[1],d0[3],0.5*(d1[4][:,0:480,:].mean(axis=1)+d1[4][:,528:1008,:].mean(axis=1)),norm=matplotlib.colors.LogNorm(),cmap="gnuplot2")
plt.colorbar()

plt.show()
