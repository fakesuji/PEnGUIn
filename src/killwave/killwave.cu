#include <unistd.h>
#include <iostream>

#include "parameters.h"
#include "structs.h"
#include "init.h"

__device__ double get_ramping_fac(double x)
{
  //f  = cpow(csin(f*hpi), 2.0);//1.0/(1.0 + 999.0*pow(1.0-f, 8.0)) - 0.001;//
  double f = __sinf(hpi*x);
  f *= f; 
  return f;
}


__global__ void killwave(Grid G, Cell* C, double dt)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;

	double rad, azi, pol;

	double f, tau;
	#ifdef shear_box
	double d_in  = sc_h*kill_width;
	double d_out = sc_h*kill_width;
	#else
	double d_in  = get_h(xmin,0.0,hpi)*kill_width;
	double d_out = get_h(xmax,0.0,hpi)*kill_width;
	#endif

	double inner = xmin+d_in;
	double outer = xmax-d_out;

	Cell C_tmp;
	Cell I_tmp;
	int ind;

	if (i>=xpad && i<G.xarr-xpad)
	if (j>=ypad && j<G.yarr-ypad)
	if (k>=zpad && k<G.zarr-zpad)
	{
		rad = G.get_xc(i);
		ind = G.get_ind(i,j,k);

		if (rad<inner) 
		{
			azi = G.get_yc(j);
			pol = G.get_zc(k);

			f = (inner-rad)/d_in;
			#if geomx == 0
			tau  = 1.0;
			#else
			tau  = pow(rad,1.5);
			#endif

			f  = get_ramping_fac(f);
			f *= dt/tau;

			C_tmp.copy(C[ind]);
			I_tmp = init_C(rad,azi,pol);

			//C_tmp.r += f * ( I_tmp.r - C_tmp.r );
			//C_tmp.p += f * ( I_tmp.p - C_tmp.p );
			C_tmp.u += f * ( I_tmp.u - C_tmp.u );
			C_tmp.v += f * ( I_tmp.v - C_tmp.v );
			C_tmp.w += f * ( I_tmp.w - C_tmp.w );

			C[ind].copy(C_tmp);
		}
		
		if (rad>outer)
		{
			azi = G.get_yc(j);
			pol = G.get_zc(k);

			f = (rad-outer)/d_out;
			#if geomx == 0
			tau  = 1.0;
			#else
			tau  = pow(rad,1.5);
			#endif

			f  = get_ramping_fac(f);
			f *= dt/tau;

			C_tmp.copy(C[ind]);
			I_tmp = init_C(rad,azi,pol);

			//C_tmp.r += f * ( I_tmp.r - C_tmp.r );
			//C_tmp.p += f * ( I_tmp.p - C_tmp.p );
			C_tmp.u += f * ( I_tmp.u - C_tmp.u );
			C_tmp.v += f * ( I_tmp.v - C_tmp.v );
			C_tmp.w += f * ( I_tmp.w - C_tmp.w );

			C[ind].copy(C_tmp);
		}
	}
	return;
}

void killwave(Grid* dev, double dt)
{
	int nx, ny, nz;
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		killwave<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xdiv,x_ydiv,x_zdiv), 0, dev[n].stream >>> (dev[n],dev[n].C,dt);
	}
	return;
}
