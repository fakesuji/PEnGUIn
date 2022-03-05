#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#include "parameters.h"
#include "structs.h"
#include "geom.h"
#include "EOS.h"
#include "init.h"
#include "orbit.h"

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

#include "kill/killwave.cu"
#include "force/forces.cu"
#include "planet/planet.cu"
#include "recon/recon.cu"
#include "boundary/boundary.cu"
#include "riemann/riemann.cu"
#include "advection/advection.cu"

__global__ void clear_flux(int mx, int my, int mz, Cell* F)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if (i<mx*my*mz) F[i].zero();
	return;
}

__global__ void sweepx(Grid G, double dt)
{
	__shared__ double xa[x_xthd+1], dx[x_xthd], xv[x_xthd];
	__shared__ double r[x_xthd*x_ydiv], p[x_xthd*x_ydiv], u[x_xthd*x_ydiv], v[x_xthd*x_ydiv], w[x_xthd*x_ydiv];

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

	r[i+x_xthd*j] = G.C[ind].r;
	p[i+x_xthd*j] = G.C[ind].p;
	u[i+x_xthd*j] = G.C[ind].u;
	v[i+x_xthd*j] = G.C[ind].v;
	w[i+x_xthd*j] = G.C[ind].w;
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

	/////////////////////////////////////////////////////

	Cell Del;
	#ifdef advec_flag
	Del = advection(geomx, xa, dx, xv, xpad, &r[x_xthd*j], &p[x_xthd*j], &u[x_xthd*j], &v[x_xthd*j], &w[x_xthd*j], 1.0, dt);
	#else
	Del =   riemann(geomx, xa, dx, xv, rad, &r[x_xthd*j], &p[x_xthd*j], &u[x_xthd*j], &v[x_xthd*j], &w[x_xthd*j], dt);
	#endif
	Del.multiply(G.get_yv(idy)*G.get_zv(idz));

	if (i>=xpad && i<x_xthd-xpad)
	{
		#if geomx == 1
		Del.v /= rad;
		#elif geomx == 2
		Del.v /= rad_cyl;
		Del.w /= rad;
		#endif
		G.F[ind].copy(Del);
	}

	return;
}


__global__ void sweepy(Grid G, double dt)
{
	__shared__ double ya[y_ythd+1], dy[y_ythd], yv[y_ythd];
	__shared__ double r[y_ythd*y_xdiv], p[y_ythd*y_xdiv], u[y_ythd*y_xdiv], v[y_ythd*y_xdiv], w[y_ythd*y_xdiv];

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

	r[i+y_ythd*j] = G.C[ind].r;
	p[i+y_ythd*j] = G.C[ind].p;
	u[i+y_ythd*j] = G.C[ind].v - G.get_rot(idx,idz);
	v[i+y_ythd*j] = G.C[ind].w;
	w[i+y_ythd*j] = G.C[ind].u;

	#if geomy == 3 || geomy == 4
	u[i+y_ythd*j] -= rad*frame_omega;
	#endif
	__syncthreads();

	/////////////////////////////////////////////////////

	Cell Del;
	Del = riemann(geomy, ya, dy, yv, rad, &r[y_ythd*j], &p[y_ythd*j], &u[y_ythd*j], &v[y_ythd*j], &w[y_ythd*j], dt);
	Del.multiply(G.get_xv(idx)*G.get_zv(idz));

	#if geomy > 2
	Del.multiply(1.0/rad);
	#endif

	if (i>=ypad && i<y_ythd-ypad)
	{
		G.F[ind].r = Del.r;
		G.F[ind].p = Del.p;
		G.F[ind].u = Del.w;
		#if geomy == 3 || geomy == 4
		G.F[ind].v = Del.u + (G.get_rot(idx,idz) + rad*frame_omega)*Del.r;
		#else
		G.F[ind].v = Del.u + G.get_rot(idx,idz)*Del.r;
		#endif
		G.F[ind].w = Del.v;
	}

	return;
}

__global__ void sweepz(Grid G, double dt)
{
	__shared__ double za[z_zthd+1], dz[z_zthd], zv[z_zthd];
	__shared__ double r[z_zthd*z_xdiv], p[z_zthd*z_xdiv], u[z_zthd*z_xdiv], v[z_zthd*z_xdiv], w[z_zthd*z_xdiv];

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
	
	r[i+z_zthd*j] = G.C[ind].r;
	p[i+z_zthd*j] = G.C[ind].p;
	u[i+z_zthd*j] = G.C[ind].w;
	v[i+z_zthd*j] = G.C[ind].u;
	w[i+z_zthd*j] = G.C[ind].v;
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
	Del = riemann(geomz, za, dz, zv, rad, &r[z_zthd*j], &p[z_zthd*j], &u[z_zthd*j], &v[z_zthd*j], &w[z_zthd*j], dt);
	Del.multiply(G.get_xv(idx)*G.get_yv(idy));

	#if geomz > 2
	Del.multiply(1.0/rad);
	#endif

	if (i>=zpad && i<z_zthd-zpad)
	{
		#if geomz == 5
		Del.w /= rad_cyl;
		#endif
		G.F[ind].r = Del.r;
		G.F[ind].p = Del.p;
		G.F[ind].u = Del.v;
		G.F[ind].v = Del.w;
		G.F[ind].w = Del.u;
	}

	return;
}


__global__ void update(Grid G, double dt)
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

		Q.copy(G.C[ind]);
		D.copy(G.F[ind]);
		D.multiply(dt);

		Q.p = get_energy(Q.r,Q.p,Q.u,Q.v,0.0);
		Q.r *= vol;
		Q.p *= Q.r;
		Q.u *= Q.r;
		Q.v *= Q.r;
		Q.w *= Q.r;

		Q.add(D);
		Q.r = fmax(Q.r,smallr);

		G.C[ind].r = Q.r/vol;
		G.C[ind].u = Q.u/Q.r;
		#if ndim>1
		G.C[ind].v = Q.v/Q.r;
		#endif
		#if ndim>2
		G.C[ind].w = Q.w/Q.r;
		#endif
		#ifndef advec_flag
		#if EOS_flag == 2
		#if internal_e_flag==0
		G.C[ind].p = fmax((Q.p-(Q.u*Q.u+Q.v*Q.v)/Q.r/2.0)*gamm/vol,smallp);
		#else
		G.C[ind].p = fmax(Q.p*gamm/vol,smallp);
		#endif
		#elif EOS_flag == 0
		G.C[ind].p = get_cs2(G.get_xc(i),G.get_yc(j),G.get_zc(k))*G.C[ind].r;
		#endif
		#endif
	}

	return;
}

__global__ void source(Grid G, double dt)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;

	double u,v,w;
	double xc,yc,zc;
	double fx,fy,fz;
	double gx,gy,gz;
	int ind;

	if (i>=xpad && i<G.xarr-xpad)
	if (j>=ypad && j<G.yarr-ypad)
	if (k>=zpad && k<G.zarr-zpad)
	{		
		ind = G.get_ind(i,j,k);
		u = G.C[ind].u;
		v = G.C[ind].v;
		w = G.C[ind].w;

		xc = G.get_xc(i);
		#if ndim > 1
		yc = G.get_yc(j);
		#else
		yc = 0.0;
		#endif
		#if ndim > 2
		zc = G.get_zc(k);
		#else
		zc = 0.0;
		#endif
		
		get_grav(xc,yc,zc,G.planets,dt,gx,gy,gz);
		get_fict(xc,yc,zc,u,v,w,fx,fy,fz);

		u += (gx+fx)*dt;
		#if ndim>1
		v += (gy+fy)*dt;
		#endif
		#if ndim>2
		w += (gz+fz)*dt;
		#endif

		G.C[ind].u = u;
		G.C[ind].v = v;
		G.C[ind].w = w;
	}

	return;
}

void solve(Grid* dev, double time, double dt)
{
	double hdt = 0.5*dt;
	int nx, ny, nz;
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		source<<< dim3((nx+1023)/512,ny,nz), dim3(1024,1,1), 0, dev[n].stream >>> (dev[n], hdt);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);

	#if ndim>2
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		boundz<<< dim3(nx,ny,2), dim3(zpad,1,1), 0, dev[n].stream >>>(dev[n]);
		sweepz<<< dim3(nz/z_zdiv,nx/z_xdiv,ny/z_ydiv), dim3(z_zthd,z_xdiv,z_ydiv), 2*sizeof(double)*z_zthd*z_xdiv*z_ydiv, dev[n].stream >>>
		      (dev[n], dt);
		update<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 0, dev[n].stream >>>
		      (dev[n], dt);
	}
	#endif

	#if ndim>1
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		boundy<<< dim3(nx,nz,2), dim3(ypad,1,1), 0, dev[n].stream >>>(dev[n]);
		sweepy<<< dim3(ny/y_ydiv,nx/y_xdiv,nz/y_zdiv), dim3(y_ythd,y_xdiv,y_zdiv), 2*sizeof(double)*y_ythd*y_xdiv*y_zdiv, dev[n].stream >>>
		      (dev[n], dt);

		update<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 0, dev[n].stream >>>
		      (dev[n], dt);
	}
	#endif

	boundx(dev);
	#ifdef OrbAdv_flag
	orb_boundx(dev);
	#endif
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		sweepx<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 2*sizeof(double)*x_xthd*x_ydiv*x_zdiv, dev[n].stream >>>
		      (dev[n], dt);

		update<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 0, dev[n].stream >>>
		      (dev[n], dt);
	}


	#if ndim>1
	#ifdef OrbAdv_flag
	set_OrbAdv(dev, dt);
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		boundy<<< dim3(nx,nz,2), dim3(ypad,1,1), 0, dev[n].stream >>>(dev[n]);
		advecy<<< dim3(ny/y_ydiv,nx/y_xdiv,nz/y_zdiv), dim3(y_ythd,y_xdiv,y_zdiv), 2*sizeof(double)*y_ythd*y_xdiv*y_zdiv, dev[n].stream >>>
		      (dev[n], dt);

		updadv<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 0, dev[n].stream >>>
		      (dev[n], dt);
	}
	#endif
	#endif

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		planet_evo<<< 1, n_planet, 0, dev[n].stream >>> (dev[n].planets, time+hdt, dt);

		source<<< dim3((nx+1023)/1024,ny,nz), dim3(1024,1,1), 0, dev[n].stream >>> (dev[n], hdt);

		#ifdef kill_flag
		killwave<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 0, dev[n].stream >>> (dev[n],dt);
		#endif
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);

	return;
}

