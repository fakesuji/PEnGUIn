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

__device__ double lim01(double a)
{
	return fmin(fmax(a,0.0),1.0);
}

__device__ void dimensionless_x(double rL, double r0, double rR, double &x, double &lx, double &ly)
{
	x = lim01((r0-rL)/(rR-rL));
	//#if recon_flag==0
	//lx = __logf(x);
	//ly = __logf(1.0-x);
	//#endif
	return;
}

__device__ double exp_lim(double x)
{
	return fmax(1.0+x,1.0e-6);
	//if (x>0.0) return 1.0+x;
	//else       return exp(x);
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
#include "dust/dust.cu"
#include "sweeps.cu"
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
	double vol;
	Cell Q;
	Cell D;
	int ind;
/*
	double fx;
	#if ndim>1
	double fy;
	#endif
	#if ndim>2
	double fz;
	#endif
*/
	if (i>=xpad && i<G.xarr-xpad)
	if (j>=ypad && j<G.yarr-ypad)
	if (k>=zpad && k<G.zarr-zpad)
	{
		ind = G.get_ind(i,j,k);
		vol = G.get_xv(i)*G.get_yv(j)*G.get_zv(k);

		Q.copy(in[ind]);
		D.copy(G.T[ind]);
		D.multiply(div*dt/vol);

/*
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
*/
		Q.p = get_energy(Q.r,Q.p,Q.u,Q.v,Q.w);
		Q.p *= Q.r;
		Q.u *= Q.r;
		Q.v *= Q.r;
		Q.w *= Q.r;

		Q.add(D);

		Q.u /= Q.r;
		Q.v /= Q.r;
		Q.w /= Q.r;
/*
		Q.u += 0.5*fx*dt;
		#if ndim>1
		Q.v += 0.5*fy*dt;
		#endif
		#if ndim>2
		Q.w += 0.5*fz*dt;
		#endif
*/
		if (Q.r<0.0)
		{
			Q.r = smallr;
			Q.p = get_cs2(G.get_xc(i),G.get_yc(j),G.get_zc(k))*Q.r;
			Q.u = in[ind].u;
			Q.v = in[ind].v;
			Q.w = in[ind].w;
			//printf("Error: negative density at %f %f %f\n",G.get_xc(i),G.get_yc(j),G.get_zc(k));
		}
		else if (Q.p<0.0)
		{
			Q.p = get_cs2(G.get_xc(i),G.get_yc(j),G.get_zc(k))*Q.r;
		}
		else
		{
			#if EOS_flag == 2
			#if internal_e_flag==0
			Q.p = Q.p*gamm - gamm*Q.r*(Q.u*Q.u+Q.v*Q.v+Q.w*Q.w)/2.0;
			if (Q.p<0.0) Q.p = get_cs2(G.get_xc(i),G.get_yc(j),G.get_zc(k))*Q.r;
			#else
			Q.p = Q.p*gamm;
			#endif
			#elif EOS_flag == 0
			Q.p = get_cs2(G.get_xc(i),G.get_yc(j),G.get_zc(k))*Q.r;
			#endif
		}

		out[ind].copy(Q);
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

	boundx(dev);
	#ifdef dust_flag
	boundx_dust(dev);
	#endif

	apply_source_terms_inplace(dev, 0.0, 0.0, hdt);

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

	#ifdef OrbAdv_flag
	set_OrbAdv(dev,dt);
	shift_OrbAdv(dev);
	advecty(dev,dt);
	#ifdef dust_flag
	advecty_dust(dev,dt);
	#endif
	#endif

	evolve_planet(dev,time+dt,dt);

	#ifdef visc_flag
	apply_source_terms(dev, 0.0, hdt, hdt);
	viscosity_tensor_evaluation2(dev);
	apply_source_terms_inplace(dev, 0.0, hdt, hdt);
	#else
	apply_source_terms_inplace(dev, 0.0, hdt, hdt);
	#endif



	return;
}

void solve(Grid* dev, double time, double dt)
{
	#ifdef RadPres_flag
	compute_extinction(dev, 1.0);
	#endif

	DS(dev,time,dt);

	#ifdef kill_flag
	killwave(dev, dt);
	#endif

	return;
}

