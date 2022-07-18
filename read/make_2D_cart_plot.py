import numpy as np
import matplotlib.pyplot as plt
import read_penguin as rd

xmax = 720
ymax = 1536
#label = 'h50_p3J_a-40_OA_St-20'
label = 'h25_p12J_a-10_OA'
#xmax = 576
#ymax = 1024

frame = 5

p0 = rd.load_2D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,label,frame)
p1 = rd.load_2D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,label,frame//2)
p2 = rd.load_2D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,label,frame//4)
p3 = rd.load_2D_dust_data('/mnt/penguin/fung/p2/',xmax,ymax,label,frame//8)

fig,axs = rd.make_figure_2D(p0,p1,p2,p3)
fig.savefig('CITau.png', bbox_inches='tight')
