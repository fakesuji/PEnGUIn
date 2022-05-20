#include <unistd.h>
#include <iostream>

#include "parameters.h"
#include "structs.h"
#include "EOS.h"
#include "init.h"
#include "orbit.h"
#include "planet.h"

__device__ int glo_idx(int i, int j, int k, int imax, int jmax)
{
	return i + imax*(j+jmax*k);
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

__global__ void sweepx(Grid G, Cell* C, double dt)
{
	__shared__ double xa[x_xthd+1], dx[x_xthd], xv[x_xthd];
	__shared__ double r[x_xthd*x_ydiv], p[x_xthd*x_ydiv], u[x_xthd*x_ydiv], v[x_xthd*x_ydiv], w[x_xthd*x_ydiv];
	double force;

	int i = threadIdx.x;
	int idx = i + blockIdx.x*x_xdiv;

	int j = threadIdx.y;
	int idy = j + blockIdx.y*x_ydiv + ypad;

	int k = threadIdx.z;
	int idz = k + blockIdx.z*x_zdiv + zpad;

	int ind = G.get_ind(idx,idy,idz);

	if (j==0)
	{
		xa[i] = G.get_xa(idx);
		if (i==blockDim.x-1) xa[i+1] = G.get_xa(idx+1);
		xv[i] = G.get_xv(idx);
	}
	__syncthreads();

	r[i+x_xthd*j] = C[ind].r;
	p[i+x_xthd*j] = C[ind].p;
	u[i+x_xthd*j] = C[ind].u;
	v[i+x_xthd*j] = C[ind].v;
	w[i+x_xthd*j] = C[ind].w;
	if (i>0) force = 0.5*(G.fx[ind] + G.fx[G.get_ind(idx-1,idy,idz)]);

	__syncthreads();

	if (j==0) dx[i] = xa[i+1] - xa[i];

	double rad = 0.5*(xa[i+1]+xa[i]);
	#if geomx == 1
	v[i+x_xthd*j] *= rad;
	#elif geomx == 2
	double rad_cyl = rad * sin(G.get_zc(idz));
	v[i+x_xthd*j] *= rad_cyl;
	w[i+x_xthd*j] *= rad;
	#endif
	__syncthreads();

//if (idy==2) printf("%f %f %f\n",rad,r[y_ythd*j+i],p[y_ythd*j+i]);

	/////////////////////////////////////////////////////
	Cell Del;
	Del =   riemann(geomx, xa, dx, xv, rad, &r[x_xthd*j], &p[x_xthd*j], &u[x_xthd*j], &v[x_xthd*j], &w[x_xthd*j], force, dt);
//if (idy==ypad && i>=xpad && i<x_xthd-xpad && idx>270 && idx<290) printf("%f %f: %f, %f, %e, %f, %f, %e\n",xa[i],xv[i+1]-xv[i],r[i+x_xthd*j],p[i+x_xthd*j],u[i+x_xthd*j],force,G.get_rot(idx,idz),Del.r);

	Del.multiply(G.get_yv(idy)*G.get_zv(idz));

	if (i>=xpad && i<x_xthd-xpad)
	{
		#if geomx == 1
		Del.v /= rad;
		#elif geomx == 2
		Del.v /= rad_cyl;
		Del.w /= rad;
		#endif
		G.F[ind].add(Del);
	}
	//if (idx==12 && G.get_j_shf(idx,idy,idz)==17) printf("sweepx: rad=%f azi=%f\n",G.get_xc(idx),G.get_yc(idy));

	return;
}


__global__ void sweepy(Grid G, Cell* C, double dt)
{
	__shared__ double ya[y_ythd+1], dy[y_ythd], yv[y_ythd];
	__shared__ double r[y_ythd*y_xdiv], p[y_ythd*y_xdiv], u[y_ythd*y_xdiv], v[y_ythd*y_xdiv], w[y_ythd*y_xdiv];
	double force;

	int i = threadIdx.x;
	int idy = i + blockIdx.x*y_ydiv;

	int j = threadIdx.y;
	int idx = j + blockIdx.y*y_xdiv + xpad;

	int k = threadIdx.z;
	int idz = k + blockIdx.z*x_zdiv + zpad;

	int ind = idx + G.xarr*idy + G.xarr*G.yarr*idz;

	if (j==0)
	{
		ya[i] = G.get_ya(idy);
		if (i==blockDim.x-1) ya[i+1] = G.get_ya(idy+1);
		yv[i] = G.get_yv(idy);
	}
	__syncthreads();

	if (j==0) dy[i] = ya[i+1] - ya[i];

	double rad;	
	#if geomy == 3
	rad = G.get_xc(idx);
	#elif geomy == 4
	rad = G.get_xc(idx);
	rad *= sin(G.get_zc(idz));
	#else
	rad = 1.0;
	#endif

	r[i+y_ythd*j] = C[ind].r;
	p[i+y_ythd*j] = C[ind].p;
	u[i+y_ythd*j] = C[ind].v - G.get_rot(idx,idz);
	v[i+y_ythd*j] = C[ind].w;
	w[i+y_ythd*j] = C[ind].u;

	#if geomy == 3 || geomy == 4
	u[i+y_ythd*j] -= rad*frame_omega;
	#endif

	if (i>0) force = 0.5*(G.fy[ind] + G.fy[G.get_ind(idx,idy-1,idz)]);
	__syncthreads();

	/////////////////////////////////////////////////////

	Cell Del;
	Del = riemann(geomy, ya, dy, yv, rad, &r[y_ythd*j], &p[y_ythd*j], &u[y_ythd*j], &v[y_ythd*j], &w[y_ythd*j], force, dt);

	Del.multiply(G.get_xv(idx)*G.get_zv(idz));

	#if geomy > 2
	Del.multiply(1.0/rad);
	#endif

	if (i>=ypad && i<y_ythd-ypad)
	{
		G.F[ind].r += Del.r;
		G.F[ind].p += Del.p;
		G.F[ind].u += Del.w;
		#if geomy == 3 || geomy == 4
		G.F[ind].v += Del.u + (G.get_rot(idx,idz) + rad*frame_omega)*Del.r;
		#else
		G.F[ind].v += Del.u + G.get_rot(idx,idz)*Del.r;
		#endif
		G.F[ind].w += Del.v;
	}

	return;
}

__global__ void sweepz(Grid G, Cell* C, double dt)
{
	__shared__ double za[z_zthd+1], dz[z_zthd], zv[z_zthd];
	__shared__ double r[z_zthd*z_xdiv], p[z_zthd*z_xdiv], u[z_zthd*z_xdiv], v[z_zthd*z_xdiv], w[z_zthd*z_xdiv];
	double force;

	int i = threadIdx.x;
	int idz = i + blockIdx.x*z_zdiv;

	int j = threadIdx.y;
	int idx = j + blockIdx.y*z_xdiv + xpad;

	int k = threadIdx.z;
	int idy = k + blockIdx.z*z_ydiv + ypad;

	int ind = G.get_ind(idx,idy,idz);

	if (j==0)
	{
		za[i] = G.get_za(idz);
		if (i==blockDim.x-1) za[i+1] = G.get_za(idz+1);
		zv[i] = G.get_zv(idz);
	}
	__syncthreads();

	if (j==0) dz[i] = za[i+1] - za[i];
	
	r[i+z_zthd*j] = C[ind].r;
	p[i+z_zthd*j] = C[ind].p;
	u[i+z_zthd*j] = C[ind].w;
	v[i+z_zthd*j] = C[ind].u;
	w[i+z_zthd*j] = C[ind].v;
	if (i>0) force = 0.5*(G.fz[ind] + G.fz[G.get_ind(idx,idy,idz-1)]);
	__syncthreads();

	double rad;
	rad = G.get_xc(idx);
	#if geomz == 5
	double rad_cyl;
	rad_cyl = rad * sin(0.5*(za[i+1]+za[i]));
	w[i+z_zthd*j] *= rad_cyl;
	#endif
	__syncthreads();

	/////////////////////////////////////////////////////

	Cell Del;
	Del = riemann(geomz, za, dz, zv, rad, &r[z_zthd*j], &p[z_zthd*j], &u[z_zthd*j], &v[z_zthd*j], &w[z_zthd*j], force, dt);
	Del.multiply(G.get_xv(idx)*G.get_yv(idy));

	#if geomz > 2
	Del.multiply(1.0/rad);
	#endif

	if (i>=zpad && i<z_zthd-zpad)
	{
		#if geomz == 5
		Del.w /= rad_cyl;
		#endif
		G.F[ind].r += Del.r;
		G.F[ind].p += Del.p;
		G.F[ind].u += Del.v;
		G.F[ind].v += Del.w;
		G.F[ind].w += Del.u;
	}

	return;
}

__global__ void update(Grid G, Cell* in, Cell* out, double dt, double div=1.0, int axis=0)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;
	double vol;
	Cell Q;
	Cell D;
	int ind;

	if (i>=xpad && i<G.xarr-xpad)
	if (j>=ypad && j<G.yarr-ypad)
	if (k>=zpad && k<G.zarr-zpad)
	{
		ind = i + G.xarr*j + G.xarr*G.yarr*k;
		vol = G.get_xv(i)*G.get_yv(j)*G.get_zv(k);

		Q.copy(in[ind]);
		D.copy(G.F[ind]);
		D.multiply(div*dt);

		Q.p = get_energy(Q.r,Q.p,Q.u,Q.v,Q.w);
		Q.r *= vol;
		Q.p *= Q.r;
		Q.u *= Q.r;
		Q.v *= Q.r;
		Q.w *= Q.r;

		if (axis==0 || axis==1) Q.u += 0.5*Q.r*G.fx[ind]*dt;
		if (axis==0 || axis==2) Q.v += 0.5*Q.r*G.fy[ind]*dt;
		if (axis==0 || axis==3) Q.w += 0.5*Q.r*G.fz[ind]*dt;

		Q.add(D);

		if (axis==0 || axis==1) Q.u += 0.5*Q.r*G.fx[ind]*dt;
		if (axis==0 || axis==2) Q.v += 0.5*Q.r*G.fy[ind]*dt;
		if (axis==0 || axis==3) Q.w += 0.5*Q.r*G.fz[ind]*dt;

		Q.u /= Q.r;
		Q.v /= Q.r;
		Q.w /= Q.r;
		Q.r /= vol;

		if (Q.r<=0.0)
		{
			Q.r = smallr;
			Q.p = get_cs2(G.get_xc(i),G.get_yc(j),G.get_zc(k))*Q.r;
			Q.u = in[ind].u;
			Q.v = in[ind].v;
			Q.w = in[ind].w;
			printf("Error: negative density at %f %f %f\n",G.get_xc(i),G.get_yc(j),G.get_zc(k));
		}
		else
		{
			#if EOS_flag == 2
			#if internal_e_flag==0
			Q.p = fmax(Q.p*gamm/vol - gamm*Q.r*(Q.u*Q.u+Q.v*Q.v+Q.w*Q.w)/2.0,smallp);
			#else
			Q.p = fmax(Q.p*gamm/vol,smallp);
			#endif
			#elif EOS_flag == 0
			Q.p = get_cs2(G.get_xc(i),G.get_yc(j),G.get_zc(k))*Q.r;
			#endif
		}

		out[ind].copy(Q);
	}

	return;
}

__global__ void sourcex(Grid G, Cell* C, double dt=0.0, double div=0.0)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;

	Cell T;
	double u,v,w;
	double xc,yc,zc;
	double fx;
	int ind;

	if (i>=1 && i<G.xarr-1)
	if (j>=ypad && j<G.yarr-ypad)
	if (k>=zpad && k<G.zarr-zpad)
	{		
		ind = G.get_ind(i,j,k);
		T.copy(C[ind]);
		u = T.u;
		v = T.v;
		w = T.w;

		xc = G.get_xc(i);
		#if ndim > 1
		yc = G.get_yc(j) + G.get_rot(i,k)*dt;
		#else
		yc = 0.0;
		#endif
		#if ndim > 2
		zc = G.get_zc(k);
		#else
		zc = 0.0;
		#endif

		fx = get_fx(xc,yc,zc,u,v,w,G.planets);
		//#ifdef visc_flag
		//fx += viscous_fx(G, C, i, j, k)/T.r;
		//#endif

		G.fx[ind] = G.fx[ind]*div + fx*(1.0-div);
	}

	return;
}

__global__ void sourcey(Grid G, Cell* C, double dt=0.0, double div=0.0)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;

	Cell T;
	double u,v,w;
	double xc,yc,zc;
	double fy;
	int ind;

	if (i>=xpad && i<G.xarr-xpad)
	if (j>=1 && j<G.yarr-1)
	if (k>=zpad && k<G.zarr-zpad)
	{		
		ind = G.get_ind(i,j,k);
		T.copy(C[ind]);	
		u = T.u;
		v = T.v;
		w = T.w;

		xc = G.get_xc(i);
		#if ndim > 1
		yc = G.get_yc(j) + G.get_rot(i,k)*dt;
		#else
		yc = 0.0;
		#endif
		#if ndim > 2
		zc = G.get_zc(k);
		#else
		zc = 0.0;
		#endif

		fy = get_fy(xc,yc,zc,u,v,w,G.planets);
		//#ifdef visc_flag
		//fy += viscous_fy(G, C, i, j, k)/T.r;
		//#endif

		G.fy[ind] = G.fy[ind]*div + fy*(1.0-div);
	}

	return;
}

__global__ void sourcez(Grid G, Cell* C, double dt=0.0, double div=0.0)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z;

	Cell T;
	double u,v,w;
	double xc,yc,zc;
	double fz;
	int ind;

	if (i>=xpad && i<G.xarr-xpad)
	if (j>=ypad && j<G.yarr-ypad)
	if (k>=1 && k<G.zarr-1)
	{		
		ind = G.get_ind(i,j,k);
		T.copy(C[ind]);	
		u = T.u;
		v = T.v;
		w = T.w;

		xc = G.get_xc(i);
		#if ndim > 1
		yc = G.get_yc(j) + G.get_rot(i,k)*dt;
		#else
		yc = 0.0;
		#endif
		#if ndim > 2
		zc = G.get_zc(k);
		#else
		zc = 0.0;
		#endif

		fz = get_fz(xc,yc,zc,u,v,w,G.planets);
		//#ifdef visc_flag
		//fz += viscous_fz(G, C, i, j, k)/T.r;
		//#endif

		G.fz[ind] = G.fz[ind]*div + fz*(1.0-div);
	}

	return;
}

__global__ void apply_viscosity(Grid G, Cell* C, double dt)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;

	Cell T;
	double fx;
	#if ndim > 1
	double fy;
	#endif
	#if ndim > 2
	double fz;
	#endif
	int ind;

	if (i>=xpad && i<G.xarr-xpad)
	if (j>=ypad && j<G.yarr-ypad)
	if (k>=zpad && k<G.zarr-zpad)
	{		
		ind = G.get_ind(i,j,k);
		T.copy(C[ind]);

		fx = viscous_fx(G, C, i, j, k)/T.r;
		T.u += fx*dt;
		#if ndim > 1
		fy = viscous_fy(G, C, i, j, k)/T.r;
		T.v += fy*dt;
		#endif
		#if ndim > 2
		fz = viscous_fz(G, C, i, j, k)/T.r;
		T.w += fz*dt;
		#endif

		C[ind].copy(T);
	}

	return;
}

__global__ void compute_forces(Grid G, Cell* C, double dt)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;

	Cell T;
	double xc,yc,zc;
	double fx;
	#if ndim > 1
	double fy;
	#endif
	#if ndim > 2
	double fz;
	#endif
	int ind;

	if (i>=xpad && i<G.xarr-xpad)
	if (j>=ypad && j<G.yarr-ypad)
	if (k>=zpad && k<G.zarr-zpad)
	{		
		ind = G.get_ind(i,j,k);
		T.copy(C[ind]);

		xc = G.get_xc(i);
		#if ndim > 1
		yc = G.get_yc(j) + G.get_rot(i,k)*0.5*dt;
		#else
		yc = 0.0;
		#endif
		#if ndim > 2
		zc = G.get_zc(k);
		#else
		zc = 0.0;
		#endif

		fx = get_fx(xc,yc,zc,T.u,T.v,T.w,G.planets);
		#ifdef visc_flag
		fx += viscous_fx(G, C, i, j, k)/T.r;
		#endif
		G.fx[ind] = fx;

		#if ndim > 1
		fy = get_fy(xc,yc,zc,T.u,T.v,T.w,G.planets);
		#ifdef visc_flag
		fy += viscous_fy(G, C, i, j, k)/T.r;
		#endif
		G.fy[ind] = fy;
		#endif

		#if ndim > 2
		fz = get_fz(xc,yc,zc,T.u,T.v,T.w,G.planets);
		#ifdef visc_flag
		fz += viscous_fz(G, C, i, j, k)/T.r;
		#endif
		G.fz[ind] = fz;
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
	double xc,yc,zc;
	double fx;
	#if ndim > 1
	double fy;
	#endif
	#if ndim > 2
	double fz;
	#endif
	int ind;

	if (i>=xpad && i<G.xarr-xpad)
	if (j>=ypad && j<G.yarr-ypad)
	if (k>=zpad && k<G.zarr-zpad)
	{		
		ind = G.get_ind(i,j,k);
		T.copy(C[ind]);

		xc = G.get_xc(i);
		#if ndim > 1
		yc = G.get_yc(j) + G.get_rot(i,k)*0.5*dt;
		#else
		yc = 0.0;
		#endif
		#if ndim > 2
		zc = G.get_zc(k);
		#else
		zc = 0.0;
		#endif

		fx = get_fx(xc,yc,zc,T.u,T.v,T.w,G.planets);
		#ifdef visc_flag
		fx += viscous_fx(G, C, i, j, k)/T.r;
		#endif
		T.u += fx*dt;

		#if ndim > 1
		fy = get_fy(xc,yc,zc,T.u,T.v,T.w,G.planets);
		#ifdef visc_flag
		fy += viscous_fy(G, C, i, j, k)/T.r;
		#endif
		T.v += fy*dt;
		#endif

		#if ndim > 2
		fz = get_fz(xc,yc,zc,T.u,T.v,T.w,G.planets);
		#ifdef visc_flag
		fz += viscous_fz(G, C, i, j, k)/T.r;
		#endif
		T.w += fz*dt;
		#endif

		C[ind].copy(T);
	}

	return;
}

void RK2(Grid* dev, double time, double dt)
{
	double hdt = 0.5*dt;
	int nx, ny, nz;
	int mx, my, mz;
	int bsz = 1024;

	syncallstreams(dev);

	boundx(dev);
	#if ndim>1
	boundy(dev);
	#endif
	#if ndim>2
	boundz(dev);
	#endif

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
		clear_forces<<< (mx*my*mz+bsz-1)/bsz, bsz, 0, dev[n].stream >>>(mx, my, mz, dev[n].fx, dev[n].fy, dev[n].fz);

		sweepx<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 2*sizeof(double)*x_xthd*x_ydiv*x_zdiv, dev[n].stream >>>
		      (dev[n], dev[n].C, hdt);

		#if ndim>1
		sweepy<<< dim3(ny/y_ydiv,nx/y_xdiv,nz/y_zdiv), dim3(y_ythd,y_xdiv,y_zdiv), 2*sizeof(double)*y_ythd*y_xdiv*y_zdiv, dev[n].stream >>>
		      (dev[n], dev[n].C, hdt);
		#endif

		#if ndim>2
		sweepz<<< dim3(nz/z_zdiv,nx/z_xdiv,ny/z_ydiv), dim3(z_zthd,z_xdiv,z_ydiv), 2*sizeof(double)*z_zthd*z_xdiv*z_ydiv, dev[n].stream >>>
		      (dev[n], dev[n].C, hdt);
		#endif

		update<<< dim3(nx/x_xdiv,ny,nz), x_xthd, 0, dev[n].stream >>> (dev[n], dev[n].C, dev[n].C, hdt);
	}

	#ifdef visc_flag
	viscosity_tensor_evaluation1(dev);
	#endif

	evolve_planet(dev,time+hdt,hdt);

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		mx = dev[n].xarr;
		my = dev[n].yarr;
		mz = dev[n].zarr;

		apply_forces<<< dim3(nx/x_xdiv,ny,nz), x_xthd, 0, dev[n].stream >>> (dev[n], dev[n].C, dt);
	}

	#ifdef kill_flag
	killwave(dev, dt);
	#endif

	#ifdef cool_flag
	cooling(dev, dt);
	#endif

	syncallstreams(dev);

	boundx(dev);
	#if ndim>1
	boundy(dev);
	#endif
	#if ndim>2
	boundz(dev);
	#endif

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
		      (dev[n], dev[n].C, hdt);

		#if ndim>1
		sweepy<<< dim3(ny/y_ydiv,nx/y_xdiv,nz/y_zdiv), dim3(y_ythd,y_xdiv,y_zdiv), 2*sizeof(double)*y_ythd*y_xdiv*y_zdiv, dev[n].stream >>>
		      (dev[n], dev[n].C, hdt);
		#endif

		#if ndim>2
		sweepz<<< dim3(nz/z_zdiv,nx/z_xdiv,ny/z_ydiv), dim3(z_zthd,z_xdiv,z_ydiv), 2*sizeof(double)*z_zthd*z_xdiv*z_ydiv, dev[n].stream >>>
		      (dev[n], dev[n].C, hdt);
		#endif

		update<<< dim3(nx/x_xdiv,ny,nz), x_xthd, 0, dev[n].stream >>> (dev[n], dev[n].C, dev[n].C, hdt);
	}

	evolve_planet(dev,time+dt,hdt);

	return;
}

void dimensional_splitting(Grid* dev, double time, double dt)
{
	double hdt = 0.5*dt;
	int nx, ny, nz;
	int mx, my, mz;
	int bsz = 1024;

	evolve_planet(dev,time+hdt,hdt);

	#ifdef kill_flag
	killwave(dev, hdt);
	#endif

	#ifdef cool_flag
	cooling(dev,hdt);
	#endif

	syncallstreams(dev);

	#ifdef visc_flag
	viscosity_tensor_evaluation1(dev);
	#else
	boundx(dev);
	#endif
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

		sourcex<<< dim3((dev[n].xarr+bsz-1)/bsz,ny,nz), bsz, 0, dev[n].stream >>> (dev[n], dev[n].C, hdt, 0.0);
		sweepx<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 2*sizeof(double)*x_xthd*x_ydiv*x_zdiv, dev[n].stream >>>
		      (dev[n], dev[n].C, dt);

		update<<< dim3(nx/x_xdiv,ny,nz), x_xthd, 0, dev[n].stream >>> (dev[n], dev[n].C, dev[n].C, dt, 1.0, 1);
	}

	#if ndim>2

	syncallstreams(dev);

	#ifdef visc_flag
	viscosity_tensor_evaluation1(dev);
	#else
	boundz(dev);
	#endif

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

		sourcez<<< dim3((nx+bsz-1)/bsz,ny,dev[n].zarr), bsz, 0, dev[n].stream >>> (dev[n], dev[n].C, hdt, 0.0);
		sweepz<<< dim3(nz/z_zdiv,nx/z_xdiv,ny/z_ydiv), dim3(z_zthd,z_xdiv,z_ydiv), 2*sizeof(double)*z_zthd*z_xdiv*z_ydiv, dev[n].stream >>>
		      (dev[n], dev[n].C, dt);

		update<<< dim3(nx/x_xdiv,ny,nz), x_xthd, 0, dev[n].stream >>> (dev[n], dev[n].C, dev[n].C, dt, 1.0, 3);
	}
	#endif

	#if ndim>1

	syncallstreams(dev);

	#ifdef visc_flag
	viscosity_tensor_evaluation1(dev);
	#else
	boundy(dev);
	#endif

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

		sourcey<<< dim3((nx+bsz-1)/bsz,dev[n].yarr,nz), bsz, 0, dev[n].stream >>> (dev[n], dev[n].C, hdt, 0.0);
		sweepy<<< dim3(ny/y_ydiv,nx/y_xdiv,nz/y_zdiv), dim3(y_ythd,y_xdiv,y_zdiv), 2*sizeof(double)*y_ythd*y_xdiv*y_zdiv, dev[n].stream >>>
		      (dev[n], dev[n].C, dt);

		update<<< dim3(nx/x_xdiv,ny,nz), x_xthd, 0, dev[n].stream >>> (dev[n], dev[n].C, dev[n].C, dt, 1.0, 2);
	}
	#endif

	evolve_planet(dev,time+dt,hdt);

	#ifdef cool_flag
	cooling(dev,hdt);
	#endif

	#ifdef kill_flag
	killwave(dev, hdt);
	#endif

	return;
}

void solve(Grid* dev, double time, double dt)
{
	double hdt = 0.5*dt;
	
	#ifndef advec_flag

		#ifdef RadPres_flag
		compute_extinction(dev, 1.0);
		#endif

		#if mode_flag == 0
		RK2(dev,time,dt);
		#elif mode_flag == 1
		dimensional_splitting(dev,time,dt);
		#endif

		#ifdef OrbAdv_flag
		set_OrbAdv(dev,dt);
		shift_OrbAdv(dev);
		advecty(dev,dt);
		#endif
	#else

		advectx(dev,dt);

		#if ndim>1
		advecty(dev,dt);
		#endif

	#endif

	return;
}

