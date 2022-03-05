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

__global__ void init_orbital_advection(double* xa, double* za, int mx, int mz, double* orb_rot, int* orb_shf)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int k = threadIdx.y + blockIdx.y*blockDim.y;
	int ind = i + mx*k;

	double rad = 0.5*(xa[i+1]+xa[i]);
	double rad_cyl = rad * sin(0.5*(za[k+1]+za[k]));

	#ifdef OrbAdv_flag
	orb_rot[ind] = get_v(rad,0.0,0.5*(za[k+1]+za[k])) - rad_cyl*frame_omega;
	#else
	orb_rot[ind] = 0.0;
	#endif
	orb_shf[ind] = 0;

	return;
}

__global__ void set_orbital_advection(Grid G, double dy, double dt)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int k = threadIdx.y + blockIdx.y*blockDim.y + zpad;

	double rad = G.get_xc(i);
	double rad_cyl = rad * sin(G.get_zc(k));
	dy *= rad_cyl;

	double rot = G.get_rot(i,k);
	int shift = __double2int_rd(rot*dt/dy + 0.5);
	double speed = shift*dy/dt;
	
	speed = rot - speed;
	shift = orb_jlim(G.get_shf(i,k)+shift,G.yres);

	G.write_res(i,k,speed);
	G.write_shf(i,k,shift);
	return;
}

void init_OrbAdv(Grid* dev)
{
	for (int i=0; i<ndev; i++) 
	{
		cudaSetDevice(i);
		init_orbital_advection<<< dim3(dev[i].xarr,dev[i].zarr,1), dim3(1,1,1), 0, dev[i].stream >>>
                                      (&dev[i].xa[dev[i].xbgn], dev[i].za, dev[i].xarr, dev[i].zarr, dev[i].orb_rot, dev[i].orb_shf);
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

