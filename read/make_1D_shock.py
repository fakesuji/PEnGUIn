import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd

label_moc = 'h50_MOC'
label_van = 'h50_VAN'
label_pe3 = 'h50_PEM3'
label_pp3 = 'h50_PPM3'
label_pe4 = 'h50_PEM4'
label_pp4 = 'h50_PPM4'

xmax = 9600
frame = 380

van = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_van,frame)
moc = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_moc,frame)
pe3 = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_pe3,frame)
pp3 = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_pp3,frame)
pe4 = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_pe4,frame)
pp4 = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_pp4,frame)

pe3_r = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_pe3+'_rev',frame)
pp3_r = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_pp3+'_rev',frame)

#pem_r = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_pem+'_rev',frame)
#plm_r = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_plm+'_rev',frame)
#ppm_r = rd.load_1D_data('/mnt/penguin/fung/p2/',xmax,label_ppm+'_rev',frame)
#a = np.flip(pem_r[2])
#b = np.flip(plm_r[2])
#c = np.flip(ppm_r[2])

xc = van[1]
plt.figure(1)


#plt.plot(xc,van[2],'go', mfc='none',label="VAN")
#plt.plot(xc,moc[2],'co', mfc='none',label="MOC")
plt.plot(xc,pe3[2],'r^', mfc='none',label="PEM3")
plt.plot(xc,pp3[2],'bs', mfc='none',label="PPM3")
#plt.plot(xc,pe4[2],'m^', mfc='none',label="PEM4")
#plt.plot(xc,pp4[2],'ks', mfc='none',label="PPM4")

plt.plot(xc,np.flip(pe3_r[2]),'r-', mfc='none',label="PEM3_rev")
plt.plot(xc,np.flip(pp3_r[2]),'b-', mfc='none',label="PPM3_rev")

print(np.std(np.flip(pe3_r[2])-pe3[2]))
print(np.std(np.flip(pp3_r[2])-pp3[2]))

plt.legend()
plt.show()
