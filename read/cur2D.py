import read as pg
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

d1 = pg.load_2D_data("/mnt/penguin/fung/p2/",864,2016,"xRange_17_250_OA_PLM",1)
m1 = pg.load_2D_data("/mnt/penguin/fung/p2/",864,2016,"xRange_17_250_OA_PLM_multi",1)
print(np.max(np.abs(d1[3]-m1[3])))
#print(d1[3][0,:]-s1[3][0,:])
#print(d1[3][0,:]-m1[3][0,:])

plt.figure(1)
plt.pcolormesh(d1[1],d1[2],d1[3],norm=matplotlib.colors.LogNorm())
#plt.pcolormesh(d1[1],d1[2],d1[3],vmax=5.0,vmin=0.0)
plt.colorbar()

plt.figure(2)
plt.pcolormesh(d1[1],d1[2],m1[3],norm=matplotlib.colors.LogNorm())
#plt.pcolormesh(d1[1],d1[2],m1[3],vmax=5.0,vmin=0.0)
plt.colorbar()

plt.figure(3)
plt.pcolormesh(d1[1],d1[2],d1[3]-m1[3])
plt.colorbar()

#plt.figure(2)
#plt.pcolormesh(d1[1],d1[2],n1[3]-d1[3])
#print(np.max(np.abs(n1[3]-d1[3])))
#print(n1[3]-d1[3])
#plt.colorbar()

#plt.figure(3)
#plt.plot(pg.cell_center(d1[1]),d1[3].mean(axis=0))
#plt.plot(pg.cell_center(d1[1]),n1[3].mean(axis=0))
#plt.yscale('log')

#plt.figure(4)
#plt.plot(pg.cell_center(d1[1]),(d1[5]*d1[3]).mean(axis=0))
#plt.plot(pg.cell_center(d0[1]),(d0[5]*d0[3]).mean(axis=0))

#plt.figure(5)
#plt.plot(pg.cell_center(d1[1]),d1[6].mean(axis=0))

plt.show()
