#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

#include "parameters.h"
#include "structs.h"
#include "util.h"

double h_get_dt(double xa, double dx, double r, double p, double u)
{
	return CFL*fabs(dx/u);
}

__global__ void get_dt_lv1(double* dt_lv1, double* xa, double* ya, double* za, int nx, int ny, int nz, Cell* C, Dust* CD, double* rot)
{
	extern __shared__ double sm[];
	int i = threadIdx.x;
	int ib = blockIdx.x;
	int ig = i + ib*blockDim.x;
	
	int nmax = (nx*ny*nz + blockDim.x*gridDim.x - 1)/(blockDim.x*gridDim.x);

	int idx, idy, idz;
	int ind;
	double cs, wid, spd;
	double tmp = 0.0;
	double rad, rad_cyl;

	for ( int n=0; n<nmax; n++)
	{
		idz = (n+nmax*ig)/(nx*ny);
		idy = (n+nmax*ig - nx*ny*idz)/(nx);
		idx = (n+nmax*ig - nx*idy - nx*ny*idz);
	
		if (idx<nx && idy<ny && idz<nz)
		{
			idx += xpad;
			idy += ypad;
			idz += zpad;
			ind = idx + (nx+2*xpad)*(idy + (ny+2*ypad)*idz);

			#if geomy == 3
			rad = 0.5*(xa[idx+1]+xa[idx]);
			rad_cyl = rad;
			#elif geomy == 4
			rad = 0.5*(xa[idx+1]+xa[idx]);
			rad_cyl = rad * sin(0.5*(za[idz+1]+za[idz]));
			#else	
			rad = 1.0;
			rad_cyl = 1.0;
			#endif

			#ifdef advec_flag
			cs = 0.0;		
			#else
			cs = sqrt(gam*C[ind].p/C[ind].r);
			#endif

			spd = fabs(C[ind].u)+cs;
			#ifdef dust_flag
			spd = fmax(fabs(CD[ind].u),spd);
			#endif
			wid = fabs(xa[idx+1]-xa[idx]);
			tmp = fmax(spd/wid,tmp);

			#if ndim>1
			spd = fabs(C[ind].v - rad_cyl*frame_omega - rot[idx+(nx+2*xpad)*idz])+cs;
			#ifdef dust_flag
			spd = fmax(fabs(CD[ind].v - rad_cyl*frame_omega - rot[idx+(nx+2*xpad)*idz]),spd);
			#endif
			wid = fabs(ya[idy+1]-ya[idy])*rad_cyl;
			tmp = fmax(spd/wid,tmp);
			#endif

			#if ndim>2
			spd = fabs(C[ind].w)+cs;
			#ifdef dust_flag
			spd = fmax(fabs(CD[ind].w),spd);
			#endif
			wid = fabs(za[idz+1]-za[idz])*rad;
			tmp = fmax(spd/wid,tmp);
			#endif
		}
	}
	
	sm[i] = tmp;
	__syncthreads();

	round_reduc_max(blockDim.x, sm);
	if (i==0) dt_lv1[ib] = sm[0];

	return;
}

__global__ void get_dt_lv2(double* dt, double* dt_lv1)
{
	extern __shared__ double sm[];
	int i = threadIdx.x;

	sm[i] = dt_lv1[i];
	__syncthreads();

	round_reduc_max(blockDim.x, sm);
	if (i==0) dt[0] = CFL/sm[0];

	return;
}

double global_dt(Grid* hst, Grid* dev, double dt)
{
	int lv1_size, nx, ny, nz;
	
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		lv1_size = min(1024,(nx*ny*nz+std_thd-1)/std_thd);
		get_dt_lv1<<< lv1_size, std_thd, std_thd*sizeof(double), dev[n].stream >>>(dev[n].Buff, &dev[n].xa[dev[n].xbgn], dev[n].ya, dev[n].za, nx, ny, nz, dev[n].C, dev[n].CD, dev[n].orb_rot);
		get_dt_lv2<<< 1, lv1_size, lv1_size*sizeof(double), dev[n].stream >>>(dev[n].dt, dev[n].Buff);
		cudaMemcpyAsync( hst[n].dt, dev[n].dt, sizeof(double), cudaMemcpyDeviceToHost, dev[n].stream );
	}
	
	for(int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);

	double tmp = dt*1.1;
	for (int n=0; n<ndev; n++)
	{
		if (*hst[n].dt<tmp) tmp = *hst[n].dt;
	}

	return tmp;
}
