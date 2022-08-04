import numpy as np
import plotly.graph_objects as go
import read_penguin as pg
from numpy import pi
from plotly.subplots import make_subplots


def get_the_slice(x,y,z, surfacecolor):
	return go.Surface(x=x, y=y, z=z, surfacecolor=surfacecolor, coloraxis='coloraxis')

def get_lims_colors(surfacecolor):# color limits for a slice
	return np.min(surfacecolor), np.max(surfacecolor)


def colorax(vmin, vmax):
	return dict(cmin=vmin,cmax=vmax)



res = 256
rB = 0.0000042/0.035/0.035
wid = 0.3
cmin = -15.0
cmax = 1.0
fig = make_subplots(rows=1, cols=2, specs=[[{"type": "scene"}, {"type": "scene"}]])
caxis=dict(colorscale='Turbo', colorbar_thickness=25, colorbar_len=0.75, **colorax(cmin, cmax))

##############################################################

imax = 448
jmax = 544
kmax = 192

time, xa, ya, za, ro, pr, u, v, w = pg.load_3D_data_alt("/home/fung/penguin2_bk/files/CPD_448x544x192_c35_p1_s0_r100_i96_o103_flatvor_iso_rot",imax,jmax,kmax)

#imax = 288
#jmax = 448
#kmax = 128

#time, xa, ya, za, ro, pr, u, v, w = pg.load_3D_data('/mnt/penguin/fung/p2/',imax,jmax,kmax,'h100_p1J_PPM',20)
ro = -w/0.035#np.log10(ro)

xc = pg.cell_center(xa)
yc = pg.cell_center(ya)
zc = pg.cell_center(za)

##############################################################

xm, ym = np.meshgrid(np.linspace(0.0,wid,res), np.linspace(0.0,wid,res), indexing='ij')
zs = np.zeros([res,res])
sz = np.zeros([res,res])
for i in range(res):
	for j in range(res):
		sz[i,j] = pg.interpolate_3D_periodicY_reflectZ(1.0+xm[i,j]*rB,pi+ym[i,j]*rB,np.pi/2.0,xc,yc,zc,ro,0.0,2.*np.pi)

fig.add_trace(go.Surface(x=xm, y=ym, z=zs, surfacecolor=sz, 
                         colorscale='Turbo', colorbar_thickness=35, colorbar_len=0.5, **colorax(cmin, cmax), colorbar_x=0.45), row=1, col=1)
print(get_lims_colors(sz))
##############################################################

xm, zm = np.meshgrid(np.linspace(0.0,wid,res), np.linspace(0.0,wid,res), indexing='ij')
ys = np.zeros([res,res])
sy = np.zeros([res,res])
for k in range(res):
	for i in range(res):
		sy[i,k] = pg.interpolate_3D_periodicY_reflectZ(1.0+xm[i,k]*rB,pi,0.5*pi+zm[i,k]*rB,xc,yc,zc,ro,0.0,2.*np.pi)

fig.add_trace(go.Surface(x=xm, y=ys, z=zm, surfacecolor=sy, 
                         colorscale='Turbo', colorbar_thickness=35, colorbar_len=0.5, **colorax(cmin, cmax), colorbar_x=0.45), row=1, col=1)
print(get_lims_colors(sy))
##############################################################

ym, zm = np.meshgrid(np.linspace(0.0,wid,res), np.linspace(0.0,wid,res), indexing='ij')
xs = np.zeros([res,res])
sx = np.zeros([res,res])
for k in range(res):
	for i in range(res):
		sx[i,k] = pg.interpolate_3D_periodicY_reflectZ(1.0,pi+ym[i,k]*rB,0.5*pi+zm[i,k]*rB,xc,yc,zc,ro,0.0,2.*np.pi)

fig.add_trace(go.Surface(x=xs, y=ym, z=zm, surfacecolor=sx, 
                         colorscale='Turbo', colorbar_thickness=35, colorbar_len=0.5, **colorax(cmin, cmax), colorbar_x=0.45), row=1, col=1)
print(get_lims_colors(sx))
##############################################################

imax = 384
jmax = 480
kmax = 128
cmin = -15.0
cmax = 1.0

time, xa, ya, za, ro, pr, u, v, w = pg.load_3D_data_alt("/home/fung/penguin2_bk/files/CPD_384x480x128_c35_p1_s0_r100_i96_o103_ad_rot",imax,jmax,kmax)

#imax = 288
#jmax = 768
#kmax = 128

#cmin = 0.0
#cmax = 7.0

#time, xa, ya, za, ro, pr, u, v, w = pg.load_3D_data('/mnt/penguin/fung/p2/',imax,jmax,kmax,'h35_p14E_PPM',20)
#rB = 0.035

ro = -w/0.035#np.log10(ro)

xc = pg.cell_center(xa)
yc = pg.cell_center(ya)
zc = pg.cell_center(za)

##############################################################

xm, ym = np.meshgrid(np.linspace(0.0,wid,res), np.linspace(0.0,wid,res), indexing='ij')
zs = np.zeros([res,res])
sz = np.zeros([res,res])
for i in range(res):
	for j in range(res):
		sz[i,j] = pg.interpolate_3D_periodicY_reflectZ(1.0+xm[i,j]*rB,pi+ym[i,j]*rB,np.pi/2.0,xc,yc,zc,ro,0.0,2.*np.pi)

fig.add_trace(go.Surface(x=xm, y=ym, z=zs, surfacecolor=sz, 
                         colorscale='Turbo', colorbar_thickness=35, colorbar_len=0.5, **colorax(cmin, cmax)),
                         row=1,col=2)
print(get_lims_colors(sz))
##############################################################

xm, zm = np.meshgrid(np.linspace(0.0,wid,res), np.linspace(0.0,wid,res), indexing='ij')
ys = np.zeros([res,res])
sy = np.zeros([res,res])
for k in range(res):
	for i in range(res):
		sy[i,k] = pg.interpolate_3D_periodicY_reflectZ(1.0+xm[i,k]*rB,pi,0.5*pi+zm[i,k]*rB,xc,yc,zc,ro,0.0,2.*np.pi)

fig.add_trace(go.Surface(x=xm, y=ys, z=zm, surfacecolor=sy, 
                         colorscale='Turbo', colorbar_thickness=35, colorbar_len=0.5, **colorax(cmin, cmax)),
                         row=1,col=2)
print(get_lims_colors(sy))
##############################################################

ym, zm = np.meshgrid(np.linspace(0.0,wid,res), np.linspace(0.0,wid,res), indexing='ij')
xs = np.zeros([res,res])
sx = np.zeros([res,res])
for k in range(res):
	for i in range(res):
		sx[i,k] = pg.interpolate_3D_periodicY_reflectZ(1.0,pi+ym[i,k]*rB,0.5*pi+zm[i,k]*rB,xc,yc,zc,ro,0.0,2.*np.pi)

fig.add_trace(go.Surface(x=xs, y=ym, z=zm, surfacecolor=sx, 
                         colorscale='Turbo', colorbar_thickness=35, colorbar_len=0.5, **colorax(cmin, cmax)),
                         row=1,col=2)
print(get_lims_colors(sx))
##############################################################


fig.update_layout( width=1400, height=700, font_size=16, margin=dict(l=0,r=0,b=0,t=0,pad=0),
                   scene1_yaxis_range=[0.0,wid], 
                   scene1_xaxis_range=[0.0,wid], 
                   scene1_zaxis_range=[0.0,wid], 
                   scene1_aspectmode='data', 
                   scene1_camera=dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0), eye=dict(x=2.0, y=2.0, z=0.75)),
                   scene2_yaxis_range=[0.0,wid], 
                   scene2_xaxis_range=[0.0,wid], 
                   scene2_zaxis_range=[0.0,wid], 
                   scene2_aspectmode='data', 
                   scene2_camera=dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0), eye=dict(x=2.0, y=2.0, z=0.75)))

fig.update_annotations(font_size=32)
fig.update_layout(scene1 = dict(xaxis_title="x/r<sub>B</sub>", yaxis_title="y/r<sub>B</sub>", zaxis_title="z/r<sub>B</sub>"))
fig.update_layout(scene2 = dict(xaxis_title="x/r<sub>B</sub>", yaxis_title="y/r<sub>B</sub>", zaxis_title="z/r<sub>B</sub>"))
fig.for_each_xaxis(lambda axis: axis.title.update(font=dict(size=24)))
fig.for_each_yaxis(lambda axis: axis.title.update(font=dict(size=24)))

fig.write_image("vz.png")
