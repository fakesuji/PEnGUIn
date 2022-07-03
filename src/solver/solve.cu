#include <unistd.h>
#include <iostream>

#include "parameters.h"
#include "structs.h"
#include "EOS.h"
#include "init.h"
#include "orbit.h"
#include "planet.h"

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

void syncallstreams(Grid* dev)
{
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
}

#include "cool/cool.cu"
#include "kill/killwave.cu"
#include "force/forces.cu"
#include "recon/recon.cu"
#include "boundary/boundary.cu"

#include "force/viscosity.cu"
#include "riemann/riemann.cu"
#include "advection/advection.cu"
#include "sweeps.cu"

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

__global__ void apply_viscous_heat(Grid G, Cell* C, double dt)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;

	double H;
	int ind;

	if (i>=xpad && i<G.xarr-xpad)
	if (j>=ypad && j<G.yarr-ypad)
	if (k>=zpad && k<G.zarr-zpad)
	{		
		ind = G.get_ind(i,j,k);

		H = viscous_heat(G, C, i, j, k);
		C[ind].p += H*gamm*dt;
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
		fy += viscous_fy(G, C, i, j, k)/T.r;
		#endif
		if (ave) G.fy[ind] = (G.fy[ind] + fy)/2.0;
		else     G.fy[ind] = fy;
		#endif

		#if ndim > 2
		fz = get_fz(xc,yc,zc,T.u,T.v,T.w,G.planets);
		#ifdef visc_flag
		fz += viscous_fz(G, C, i, j, k)/T.r;
		#endif
		if (ave) G.fz[ind] = (G.fz[ind] + fz)/2.0;
		else     G.fz[ind] = fz;
		#endif

		fx = get_fx(xc,yc,zc,T.u,T.v,T.w,G.planets);
		#ifdef visc_flag
		fx += viscous_fx(G, C, i, j, k)/T.r;
		#endif
		if (ave) G.fx[ind] = (G.fx[ind] + fx)/2.0;
		else     G.fx[ind] = fx;
		//G.Du[ind] = -get_slope(geomx, G.get_xc(i-1), G.get_xc(i), G.get_xc(i+1), G.get_xc(i+2),
                //                       C[G.get_ind(i-2,j,k)].p, C[G.get_ind(i-1,j,k)].p, T.p, C[G.get_ind(i+1,j,k)].p, C[G.get_ind(i+2,j,k)].p);
		#if EOS_flag==2
		#if internal_e_flag==0
		//G.De[ind] = -get_slope(geomx, G.get_xc(i-1), G.get_xc(i), G.get_xc(i+1), G.get_xc(i+2),
                //                       C[G.get_ind(i-2,j,k)].p*C[G.get_ind(i-2,j,k)].u, C[G.get_ind(i-1,j,k)].p*C[G.get_ind(i-1,j,k)].u, T.p*T.u,
                //                       C[G.get_ind(i+1,j,k)].p*C[G.get_ind(i+1,j,k)].u, C[G.get_ind(i+2,j,k)].p*C[G.get_ind(i+2,j,k)].u);
		#elif internal_e_flag==1
		//G.De[ind] = -T.p*get_slope(geomx, G.get_xc(i-1), G.get_xc(i), G.get_xc(i+1), G.get_xc(i+2),
                //                           C[G.get_ind(i-2,j,k)].u, C[G.get_ind(i-1,j,k)].u, T.u, C[G.get_ind(i+1,j,k)].u, C[G.get_ind(i+2,j,k)].u);
  		#endif
		#endif
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
	int nx, ny, nz;
	int mx, my, mz;
	int bsz = 32;
	double hdt = 0.5*dt;

	//////////////////////////////////////////////////////////////

	#ifdef visc_flag
	syncallstreams(dev);
	viscosity_tensor_evaluation1(dev);

	#else
	syncallstreams(dev);

	boundx(dev);
	#if ndim>1
	boundy(dev,time);
	#endif
	#if ndim>2
	boundz(dev);
	#endif
	#endif

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		mx = dev[n].xarr;
		my = dev[n].yarr;
		mz = dev[n].zarr;

		compute_forces<<< dim3((mx+bsz-1)/bsz,my,mz), bsz, 0, dev[n].stream >>> (dev[n], dev[n].C, 0.0);
	}

	#ifdef cool_flag
	cooling(dev, hdt);
	#endif

	//////////////////////////////////////////////////////////////

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		mx = dev[n].xarr;
		my = dev[n].yarr;
		mz = dev[n].zarr;

		clear_flux<<< (mx*my*mz+bsz-1)/bsz, bsz, 0, dev[n].stream >>>(mx, my, mz, dev[n].F);

		sweepx<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 2*sizeof(double)*x_xthd*x_ydiv*x_zdiv, dev[n].stream >>>
		      (dev[n], dev[n].C, dt);

		#if ndim>1
		sweepy<<< dim3(ny/y_ydiv,nx/y_xdiv,nz/y_zdiv), dim3(y_ythd,y_xdiv,y_zdiv), 2*sizeof(double)*y_ythd*y_xdiv*y_zdiv, dev[n].stream >>>
		      (dev[n], dev[n].C, dt);
		#endif

		#if ndim>2
		sweepz<<< dim3(nz/z_zdiv,nx/z_xdiv,ny/z_ydiv), dim3(z_zthd,z_xdiv,z_ydiv), 2*sizeof(double)*z_zthd*z_xdiv*z_ydiv, dev[n].stream >>>
		      (dev[n], dev[n].C, dt);
		#endif
	}

	//////////////////////////////////////////////////////////////

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		update<<< dim3(nx/x_xdiv,ny,nz), x_xdiv, 0, dev[n].stream >>> (dev[n], dev[n].C, dev[n].T, hdt);
	}

	//////////////////////////////////////////////////////////////

	evolve_planet(dev,time+hdt,hdt);

	//////////////////////////////////////////////////////////////

	#ifdef visc_flag
	syncallstreams(dev);
	viscosity_tensor_evaluation2(dev);

	#else

	syncallstreams(dev);

	boundx2(dev);
	#if ndim>1
	boundy2(dev,time+hdt);
	#endif
	#if ndim>2
	boundz2(dev);
	#endif
	#endif

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		mx = dev[n].xarr;
		my = dev[n].yarr;
		mz = dev[n].zarr;

		compute_forces<<< dim3((mx+bsz-1)/bsz,my,mz), bsz, 0, dev[n].stream >>> (dev[n], dev[n].T, hdt);
	}

	//////////////////////////////////////////////////////////////

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		update<<< dim3(nx/x_xdiv,ny,nz), x_xdiv, 0, dev[n].stream >>> (dev[n], dev[n].C, dev[n].C, dt);
	}

	//////////////////////////////////////////////////////////////

	evolve_planet(dev,time+dt,hdt);

	//////////////////////////////////////////////////////////////

	#ifdef cool_flag
	cooling(dev, hdt);
	#endif

	#ifdef kill_flag
	killwave(dev, dt);
	#endif

	return;
}

void DS(Grid* dev, double time, double dt)
{
	double hdt = 0.5*dt;

	#ifdef visc_flag
	viscosity_tensor_evaluation1(dev);
	#endif

	source_terms_inplace(dev, 0.0, hdt);

	boundx(dev);
	sweepx_inplace(dev,dt);

	#if ndim>1
	boundy(dev,time+dt);
	sweepy_inplace(dev,dt);
	#endif

	evolve_planet(dev,time+dt,dt);

	source_terms_replace(dev, dt, hdt);
	for (int n=0; n<ndev; n++) dev[n].CT_change();

	#ifdef visc_flag
	viscosity_tensor_evaluation1(dev);
	#endif

	source_terms_update(dev, dt, hdt);
	for (int n=0; n<ndev; n++) dev[n].CT_change();

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
	#endif

	return;
}

