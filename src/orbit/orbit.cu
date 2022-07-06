#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "parameters.h"
#include "structs.h"
#include "init.h"

__device__ int orb_jlim(int j, int jmax)
{
	j = j%jmax;
	return j;
}

__global__ void shift_orbital_advection(Grid G)
{
	int i = blockIdx.x + xpad;
	int j;
	int k = blockIdx.y + zpad;

	int ind, shf;
	
	int jmax = G.yres;
	int len = blockDim.x;
	int loop = (jmax+len-1)/len;
	int inc = G.get_inc(i,k);

	for (int l=0; l<loop; l++)
	{
		j = threadIdx.x + l*len;
		if (j<jmax)
		{
			ind = i + G.xarr*(j+ypad) + G.xarr*G.yarr*k;
			j = (j + inc)%jmax;
			if (j<0) j += jmax;
			shf = i + G.xarr*(j+ypad) + G.xarr*G.yarr*k;
		
			G.T[shf].copy(G.C[ind]);

			#ifdef dust_flag
			G.TD[shf].copy(G.CD[ind]);
			#endif
		}
	}

	return;
}

__global__ void init_orbital_advection(Grid G)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int k = threadIdx.y + blockIdx.y*blockDim.y;
	int ind = i + G.xarr*k;

	double rad = G.get_xc(i);
	#if geomx==2
	double rad_cyl = rad * sin(G.get_zc(k));
	#else
	double rad_cyl = rad;
	#endif

	#ifdef OrbAdv_flag
	G.orb_rot[ind] = get_v(rad,0.0,G.get_zc(k)) - rad_cyl*frame_omega;
	#else
	G.orb_rot[ind] = 0.0;
	#endif
	G.orb_shf[ind] = 0;

	return;
}

__global__ void set_orbital_advection(Grid G, double dy, double dt)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int k = threadIdx.y + blockIdx.y*blockDim.y + zpad;

	double rad = G.get_xc(i);
	#if geomx==2
	double rad_cyl = rad * sin(G.get_zc(k));
	#else
	double rad_cyl = rad;
	#endif
	dy *= rad_cyl;

	double rot = G.get_rot(i,k);
	int shift = __double2int_rd(rot*dt/dy + 0.5);
	double speed = shift*dy/dt;
	
	speed = rot - speed;
	shift = shift;
	
	//printf("%i %i %f: %f %f %i\n",i,k,rad,rot,speed,shift);

	G.write_res(i,k,speed);
	G.write_inc(i,k,shift);
	return;
}

__global__ void update_orbital_advection(Grid G)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int k = threadIdx.y + blockIdx.y*blockDim.y + zpad;

	int shift = orb_jlim(G.get_shf(i,k)+G.get_inc(i,k),G.yres);
	G.write_shf(i,k,shift);
	return;
}

void init_OrbAdv(Grid* dev)
{
	for (int i=0; i<ndev; i++) 
	{
		cudaSetDevice(i);
		init_orbital_advection<<< dim3(dev[i].xarr,dev[i].zarr,1), dim3(1,1,1), 0, dev[i].stream >>> (dev[i]);
	}
	for(int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
	return;
}

void set_OrbAdv(Grid* dev, double dt)
{
	for (int i=0; i<ndev; i++) 
	{
		cudaSetDevice(i);
		set_orbital_advection<<< dim3(dev[i].xres/x_xdiv,dev[i].zres,1), dim3(x_xdiv,1,1), 0, dev[i].stream >>> (dev[i], (ymax-ymin)/(double)yres, dt);
	}
	for(int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
	return;
}

void update_OrbAdv(Grid* dev)
{
	for (int i=0; i<ndev; i++) 
	{
		cudaSetDevice(i);
		update_orbital_advection<<< dim3(dev[i].xres/x_xdiv,dev[i].zres,1), dim3(x_xdiv,1,1), 0, dev[i].stream >>> (dev[i]);
	}
	for(int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
	return;
}

void shift_OrbAdv(Grid* dev)
{
	for (int i=0; i<ndev; i++) 
	{
		cudaSetDevice(i);
		shift_orbital_advection<<< dim3(dev[i].xres,dev[i].zres,1), dim3(1024,1,1), 0, dev[i].stream >>> (dev[i]);
		dev[i].CT_change();
		#ifdef dust_flag
		dev[i].CT_D_change();
		#endif
	}
	for(int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
	return;
}

