import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd

label_moc = 'h50_MOC'
label_van = 'h50_VAN'
label_pe3 = 'h50_PEM3'
label_pp3 = 'h50_PPM3'
label_pe4 = 'h50_PEM4'
label_pp4 = 'h50_PPM4'

xmax = 128
frame = 10

van = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_van,frame)
moc = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_moc,frame)
pe3 = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_pe3,frame)
pp3 = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_pp3,frame)
pe4 = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_pe4,frame)
pp4 = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_pp4,frame)

xc = van[1]
sol = 1.0 + 1.e-6*np.sin(2.0*np.pi*xc)

plt.figure(1)

plt.plot(xc,van[2]-sol,'-', mfc='none',label='VAN')
print('VAN',np.std(van[2]-sol))
plt.plot(xc,moc[2]-sol,'-', mfc='none',label='MOC')
print('MOC',np.std(moc[2]-sol))
plt.plot(xc,pe3[2]-sol,'-', mfc='none',label='PEM3')
print('PEM3',np.std(pe3[2]-sol))
plt.plot(xc,pp3[2]-sol,'-', mfc='none',label='PPM3')
print('PPM3',np.std(pp3[2]-sol))
plt.plot(xc,pe4[2]-sol,'-', mfc='none',label='PEM4')
print('PEM4',np.std(pe4[2]-sol))
plt.plot(xc,pp4[2]-sol,'-', mfc='none',label='PPM4')
print('PPM4',np.std(pp4[2]-sol))

plt.legend()

plt.show()
