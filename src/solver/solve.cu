#include <unistd.h>
#include <iostream>

#include "parameters.h"
#include "structs.h"
#include "EOS.h"
#include "init.h"
#include "boundary.h"
#include "killwave.h"
#include "orbit.h"
#include "planet.h"
#include "viscosity.h"

__device__ double exp_lim(double x)
{
	return fmax(1.0+x,1.0e-6);
}

__device__ int jlim(int j, int jmax)
{
	j = j%jmax;
	if (j<0) j += jmax;
	return j;
}

__device__ double get_dv_dr_dev(int geom, double ra, double dr)
{
	if 	(geom==5)
	{
		if (dr<1.0e-4) return sin(ra+0.5*dr);
		else           return (cos(ra)-cos(ra+dr))/dr;
	}
	//else if (geom==4) return 1.0;	
	//else if (geom==3) return 1.0;
	else if (geom==2) return ra*ra + ra*dr + dr*dr/3.0;
	else if (geom==1) return ra + dr/2.0;
	else return 1.0;
}

__device__ double get_volume_dev(int geom, double ra, double dr)
{
	return dr*get_dv_dr_dev(geom,ra,dr);
}

__device__ double get_rad(Grid G, int i, int j, int k)
{
	#if geomx == 1
	return G.get_xc(i);
	#elif geomx == 2
	return G.get_xc(i);
	#else
	return 1.0;
	#endif
}

__device__ double get_rad_cyl(Grid G, int i, int j, int k)
{
	#if geomx == 1
	return G.get_xc(i);
	#elif geomx == 2
	return G.get_xc(i) * sin(G.get_zc(k));
	#else
	return 1.0;
	#endif
}

#include "cool/cool.cu"
#include "force/forces.cu"
#include "recon/recon.cu"

#include "riemann/riemann.cu"
#include "sweeps.cu"
#include "dust/dust.cu"
#include "advection/advection.cu"

__global__ void clear_flux(int mx, int my, int mz, Cell* F)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if (i<mx*my*mz) F[i].zero();
	return;
}

__global__ void clear_forces(int mx, int my, int mz, double* fx, double* fy, double* fz)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if (i<mx*my*mz) 
	{
		fx[i] = 0.0;
		fy[i] = 0.0;
		fz[i] = 0.0;
	}
	return;
}

__global__ void update(Grid G, Cell* in, Cell* out, double dt, double div=1.0)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;
	double vol, r_min, p_min;
	Cell Q;
	Cell D;
	int ind;

	double fx;
	#if ndim>1
	double fy;
	#endif
	#if ndim>2
	double fz;
	#endif

	if (i>=xpad && i<G.xarr-xpad)
	if (j>=ypad && j<G.yarr-ypad)
	if (k>=zpad && k<G.zarr-zpad)
	{
		ind = G.get_ind(i,j,k);
		vol = G.get_xv(i)*G.get_yv(j)*G.get_zv(k);

		Q.copy(in[ind]);
		D.copy(G.F[ind]);
		D.multiply(div*dt/vol);

		r_min = fmax(Q.r*1.0e-10,smallr);
		p_min = fmax(Q.p*1.0e-10,smallp);

		fx = G.fx[ind];
		Q.u += 0.5*fx*dt;
		#if ndim>1
		fy = G.fy[ind];
		Q.v += 0.5*fy*dt;
		#endif
		#if ndim>2
		fz = G.fz[ind];
		Q.w += 0.5*fz*dt;
		#endif

		Q.p = get_energy(Q.r,Q.p,Q.u,Q.v,Q.w);
		Q.p *= Q.r;
		Q.u *= Q.r;
		Q.v *= Q.r;
		Q.w *= Q.r;

		Q.add(D);

		Q.u /= Q.r;
		Q.v /= Q.r;
		Q.w /= Q.r;

		Q.u += 0.5*fx*dt;
		#if ndim>1
		Q.v += 0.5*fy*dt;
		#endif
		#if ndim>2
		Q.w += 0.5*fz*dt;
		#endif

		if (Q.r<r_min)
		{
			Q.r = r_min;
			Q.p = p_min;
			Q.u = in[ind].u;
			Q.v = in[ind].v;
			Q.w = in[ind].w;
			//printf("Error: negative density at %f %f %f\n",G.get_xc(i),G.get_yc(j),G.get_zc(k));
		}
		else
		{
			#if EOS_flag == 2
			#if internal_e_flag==0
			Q.p = fmax(Q.p*gamm - gamm*Q.r*(Q.u*Q.u+Q.v*Q.v+Q.w*Q.w)/2.0,p_min);
			#else
			Q.p = fmax(Q.p*gamm,p_min);
			#endif
			#elif EOS_flag == 0
			Q.p = get_cs2(G.get_xc(i),G.get_yc(j),G.get_zc(k))*Q.r;
			#endif
		}

		out[ind].copy(Q);
	}

	return;
}

__global__ void compute_forces(Grid G, Cell* C, double x_dt, bool ave=false)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;

	Cell T;
	double xc,yc,zc;
	double fx = 0.0;
	double fy = 0.0;
	double fz = 0.0;
	int ind;

	if (i>=xpad-1 && i<G.xarr-xpad+1)
	if (j>=ypad-1 && j<G.yarr-ypad+1)
	if (k>=zpad-1 && k<G.zarr-zpad+1)
	{		
		ind = G.get_ind(i,j,k);
		T.copy(C[ind]);

		xc = G.get_xc(i);
		#if ndim > 1
		yc = G.get_yc(j) + G.get_rot(i,k)*x_dt;
		#else
		yc = 0.0;
		#endif
		#if ndim > 2
		zc = G.get_zc(k);
		#else
		zc = 0.0;
		#endif

		#if ndim > 1
		fy = get_fy(xc,yc,zc,T.u,T.v,T.w,G.planets);
		#ifdef visc_flag
		fy += G.vis_tensor[ndim*ind+1]/T.r;
		#endif
		if (ave) G.fy[ind] = (G.fy[ind] + fy)/2.0;
		else     G.fy[ind] = fy;
		#endif

		#if ndim > 2
		fz = get_fz(xc,yc,zc,T.u,T.v,T.w,G.planets);
		#ifdef visc_flag
		fz += G.vis_tensor[ndim*ind+2]/T.r;
		#endif
		if (ave) G.fz[ind] = (G.fz[ind] + fz)/2.0;
		else     G.fz[ind] = fz;
		#endif

		fx = get_fx(xc,yc,zc,T.u,T.v,T.w,G.planets);
		#ifdef visc_flag
		fx += G.vis_tensor[ndim*ind+0]/T.r;
		#endif
		if (ave) G.fx[ind] = (G.fx[ind] + fx)/2.0;
		else     G.fx[ind] = fx;
	}

	return;
}

__global__ void apply_forces(Grid G, Cell* C, double dt)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;

	Cell T;
	int ind;

	if (i>=xpad && i<G.xarr-xpad)
	if (j>=ypad && j<G.yarr-ypad)
	if (k>=zpad && k<G.zarr-zpad)
	{		
		ind = G.get_ind(i,j,k);
		T.copy(C[ind]);

		T.u += G.fx[ind]*dt;

		#if ndim > 1
		T.v += G.fy[ind]*dt;
		#endif

		#if ndim > 2
		T.w += G.fz[ind]*dt;
		#endif

		C[ind].copy(T);
	}

	return;
}

void RK2(Grid* dev, double time, double dt)
{

	return;
}

void DS(Grid* dev, double time, double dt)
{
	double hdt = 0.5*dt;

	#ifdef visc_flag
	viscosity_tensor_evaluation1(dev);
	#endif

	source_terms_replace(dev, 0.0, hdt);
	#ifdef dust_flag
	source_terms_dust_replace(dev, 0.0, hdt);
	#endif
	for (int n=0; n<ndev; n++) dev[n].CT_change();
	#ifdef dust_flag
	for (int n=0; n<ndev; n++) dev[n].CT_D_change();
	#endif

	boundx(dev);
	sweepx_inplace(dev,dt);

	#if ndim>1
	boundy(dev,time+dt);
	sweepy_inplace(dev,dt);
	#endif

	#if ndim>2
	boundz(dev);
	sweepz_inplace(dev,dt);
	#endif

	#ifdef dust_flag
	boundx_dust(dev);
	sweepx_dust_inplace(dev,dt);

	#if ndim>1
	boundy_dust(dev,time+dt);
	sweepy_dust_inplace(dev,dt);
	#endif

	#if ndim>2
	boundz_dust(dev);
	sweepz_dust_inplace(dev,dt);
	#endif
	#endif

	evolve_planet(dev,time+dt,dt);

	source_terms_replace(dev, dt, hdt);
	for (int n=0; n<ndev; n++) dev[n].CT_change();
	#ifdef dust_flag
	source_terms_dust_replace(dev, dt, hdt);
	for (int n=0; n<ndev; n++) dev[n].CT_D_change();
	#endif

	#ifdef visc_flag
	viscosity_tensor_evaluation1(dev);
	#endif

	source_terms_update(dev, dt, hdt);
	for (int n=0; n<ndev; n++) dev[n].CT_change();
	#ifdef dust_flag
	source_terms_dust_update(dev, dt, hdt);
	for (int n=0; n<ndev; n++) dev[n].CT_D_change();
	#endif

	return;
}

void solve(Grid* dev, double time, double dt)
{
	#ifdef RadPres_flag
	compute_extinction(dev, 1.0);
	#endif

	#if mode_flag == 0
	DS(dev,time,dt);
	#elif mode_flag == 1
	#endif

	#ifdef OrbAdv_flag
	set_OrbAdv(dev,dt);
	shift_OrbAdv(dev);
	advecty(dev,dt);
	#ifdef dust_flag
	advecty_dust(dev,dt);
	#endif
	#endif

	#ifdef kill_flag
	killwave(dev, dt);
	#endif

	return;
}

