import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd

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

print(np.sqrt(np.sum((ppm[2]-sol)*(ppm[2]-sol))/xmax))

plt.plot(xc,plm[2]-sol,'g:', mfc='none')
#plt.plot(hpl[1],hpl[2],'gx-')
#plt.plot(xc,b,'gx')
#print(np.max(np.abs(plm[2]-b)))

print(np.sqrt(np.sum((plm[2]-sol)*(plm[2]-sol))/xmax))

plt.plot(xc,pem[2]-sol,'r:', mfc='none')
#plt.plot(hpe[1],hpe[2],'rx-')
#plt.plot(xc,a,'rx')
#print(np.max(np.abs(pem[2]-a)))

print(np.sqrt(np.sum((pem[2]-sol)*(pem[2]-sol))/xmax))

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
