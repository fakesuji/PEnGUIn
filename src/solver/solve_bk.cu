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
	__shared__ double ya[x_ydiv+1], yv[x_ydiv];
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
		xa[i] = G.xa[idx+G.xbgn];
		if (i==blockDim.x-1) xa[i+1] = G.xa[idx+1+G.xbgn];
		xv[i] = G.xv[idx+G.xbgn];
	}
	if (i==0)
	{
		ya[j] = G.ya[idy+G.ybgn];
		if (j==blockDim.y-1) ya[j+1] = G.ya[idy+1+G.ybgn];
		yv[j] = G.yv[idy+G.ybgn];
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
	double rad_cyl = rad * sin(0.5*(G.za[idz+1+G.zbgn]+G.za[idz+G.zbgn]));
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
	Del.multiply(yv[j]*G.zv[idz+G.zbgn]);

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


__global__ void sweepy(double* d_ya, double* d_yv, double* d_xa, double* d_xv, double* za, double* zv,
                       int mx, int my, Cell* C, Cell* F, double* rot, double dt)
{
	__shared__ double ya[y_ythd+1], dy[y_ythd], yv[y_ythd];
	__shared__ double xa[y_xdiv+1], xv[y_xdiv];
	__shared__ double r[y_ythd*y_xdiv], p[y_ythd*y_xdiv], u[y_ythd*y_xdiv], v[y_ythd*y_xdiv], w[y_ythd*y_xdiv];
	#if advection_flag == 1
	__shared__ double s[y_ythd*y_xdiv];
	#endif

	int i = threadIdx.x;
	int idy = i + blockIdx.x*y_ydiv;

	int j = threadIdx.y;
	int idx = j + blockIdx.y*y_xdiv + xpad;

	int k = threadIdx.z;
	int idz = k + blockIdx.z*x_zdiv + zpad;

	int ind = idx + mx*idy + mx*my*idz;

	if (j==0)
	{
		ya[i] = d_ya[idy];
		if (i==blockDim.x-1) ya[i+1] = d_ya[idy+1];
		yv[i] = d_yv[idy];
	}
	if (i==0)
	{
		xa[j] = d_xa[idx];
		if (j==blockDim.y-1) xa[j+1] = d_xa[idx+1];
		xv[j] = d_xv[idx];
	}
	__syncthreads();

	if (j==0) dy[i] = ya[i+1] - ya[i];

	double rad;	
	#if geomy == 3
	rad = 0.5*(xa[j+1]+xa[j]);
	#elif geomy == 4
	rad = 0.5*(xa[j+1]+xa[j]);
	rad *= sin(0.5*(za[idz+1]+za[idz]));
	#else
	rad = 1.0;
	#endif

	r[i+y_ythd*j] = C[ind].r;
	p[i+y_ythd*j] = C[ind].p;
	u[i+y_ythd*j] = C[ind].v - rot[idx+mx*idz];
	v[i+y_ythd*j] = C[ind].w;
	w[i+y_ythd*j] = C[ind].u;

	#if geomy == 3 || geomy == 4
	u[i+y_ythd*j] -= rad*frame_omega;
	#endif
	__syncthreads();

	/////////////////////////////////////////////////////

	Cell Del;
	Del = riemann(geomy, ya, dy, yv, rad, &r[y_ythd*j], &p[y_ythd*j], &u[y_ythd*j], &v[y_ythd*j], &w[y_ythd*j], dt);
	Del.multiply(xv[j]*zv[idz]);

	#if geomy > 2
	Del.multiply(1.0/rad);
	#endif

	if (i>=ypad && i<y_ythd-ypad)
	{
		F[ind].r = Del.r;
		F[ind].p = Del.p;
		F[ind].u = Del.w;
		#if geomy == 3 || geomy == 4
		F[ind].v = Del.u + (rot[idx+mx*idz] + rad*frame_omega)*Del.r;
		#else
		F[ind].v = Del.u + rot[idx+mx*idz]*Del.r;
		#endif
		F[ind].w = Del.v;
	}

	return;
}

__global__ void sweepz(double* d_za, double* d_zv, double* d_xa, double* d_xv, double* ya, double* yv,
                       int mx, int my, Cell* C, Cell* F, int* shift, double dt)
{
	__shared__ double za[z_zthd+1], dz[z_zthd], zv[z_zthd];
	__shared__ double xa[z_xdiv+1], xv[z_xdiv];
	__shared__ double r[z_zthd*z_xdiv], p[z_zthd*z_xdiv], u[z_zthd*z_xdiv], v[z_zthd*z_xdiv], w[z_zthd*z_xdiv];
	#if advection_flag == 1
	__shared__ double s[z_zthd*z_xdiv];
	#endif

	int i = threadIdx.x;
	int idz = i + blockIdx.x*z_zdiv;

	int j = threadIdx.y;
	int idx = j + blockIdx.y*z_xdiv + xpad;

	int k = threadIdx.z;
	int idy = k + blockIdx.z*z_ydiv + ypad;

	int idy_shf = jlim(k + blockIdx.z*z_ydiv - shift[idx+mx*idz], yres) + ypad;

	int ind = idx + mx*idy_shf + mx*my*idz;

	if (j==0)
	{
		za[i] = d_za[idz];
		if (i==blockDim.x-1) za[i+1] = d_za[idz+1];
		zv[i] = d_zv[idz];
	}
	if (i==0)
	{
		xa[j] = d_xa[idx];
		if (j==blockDim.y-1) xa[j+1] = d_xa[idx+1];
		xv[j] = d_xv[idx];
	}
	__syncthreads();

	r[i+z_zthd*j] = C[ind].r;
	p[i+z_zthd*j] = C[ind].p;
	u[i+z_zthd*j] = C[ind].w;
	v[i+z_zthd*j] = C[ind].u;
	w[i+z_zthd*j] = C[ind].v;
	__syncthreads();

	if (j==0) dz[i] = za[i+1] - za[i];

	double rad;
	rad = 0.5*(xa[j+1]+xa[j]);
	#if geomz == 5
	double rad_cyl;
	rad_cyl = rad * sin(0.5*(za[i+1]+za[i]));
	w[i+z_zthd*j] *= rad_cyl;
	#endif
	__syncthreads();

	/////////////////////////////////////////////////////

	Cell Del;
	Del = riemann(geomz, za, dz, zv, rad, &r[z_zthd*j], &p[z_zthd*j], &u[z_zthd*j], &v[z_zthd*j], &w[z_zthd*j], dt);
	Del.multiply(xv[j]*yv[idy]);

	#if geomz > 2
	Del.multiply(1.0/rad);
	#endif

	if (i>=zpad && i<z_zthd-zpad)
	{
		#if geomz == 5
		Del.w /= rad_cyl;
		#endif
		F[ind].r = Del.r;
		F[ind].p = Del.p;
		F[ind].u = Del.v;
		F[ind].v = Del.w;
		F[ind].w = Del.u;
	}

	return;
}


__global__ void update(double dt, double* xa, double* ya, double* za, double* xv, double* yv, double* zv, 
                       int mx, int my, int mz, Cell* C, Cell* F)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;
	double vol = xv[i]*yv[j]*zv[k];
	Cell Q;
	//Cell T;

	int ind;
	if (i>=xpad && i<mx-xpad)
	if (j>=ypad && j<my-ypad)
	if (k>=zpad && k<mz-zpad)
	{
		ind = i + mx*(j + my*k);

		//T.copy(C[ind]);
		Q.copy(C[ind]);
		Q.p = get_energy(Q.r,Q.p,Q.u,Q.v,0.0);
		Q.r *= vol;
		Q.p *= Q.r;
		Q.u *= Q.r;
		Q.v *= Q.r;
		Q.w *= Q.r;

		F[ind].multiply(dt);
		Q.add(F[ind]);
		Q.r = fmax(Q.r,smallr);

		C[ind].r = Q.r/vol;
		C[ind].u = Q.u/Q.r;
		#if ndim>1
		C[ind].v = Q.v/Q.r;
		#endif
		#if ndim>2
		C[ind].w = Q.w/Q.r;
		#endif
		#ifndef advec_flag
		#if EOS_flag == 2
		#if internal_e_flag==0
		C[ind].p = fmax((Q.p-(Q.u*Q.u+Q.v*Q.v)/Q.r/2.0)*gamm/vol,smallp);
		#else
		C[ind].p = fmax(Q.p*gamm/vol,smallp);
		#endif
		#elif EOS_flag == 0
		C[ind].p = get_cs2(0.5*(xa[i]+xa[i+1]),0.5*(ya[j]+ya[j+1]),0.5*(za[k]+za[k+1]))*C[ind].r;
		#endif
		#endif

		//if (j==2) printf("%f, %f, %.10e, %.10e, %.10e, %.10e\n",dt,0.5*(xa[i]+xa[i+1]), T.r-C[ind].r, T.p-C[ind].p, T.u-C[ind].u, T.v-C[ind].v);
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

		xc = 0.5*(G.xa[i+G.xbgn]+G.xa[i+1+G.xbgn]);
		#if ndim > 1
		yc = 0.5*(G.ya[j+G.ybgn]+G.ya[j+1+G.ybgn]);
		#else
		yc = 0.0;
		#endif
		#if ndim > 2
		zc = 0.5*(G.za[k+G.zbgn]+G.za[k+1+G.zbgn]);
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
	int nx, ny, nz, mx, my, mz;
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		mx = dev[n].xarr;
		my = dev[n].yarr;
		mz = dev[n].zarr;

		source<<< dim3((nx+1023)/512,ny,nz), dim3(1024,1,1), 0, dev[n].stream >>>
		      (dev[n], hdt);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);

	#if ndim>2
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		mx = dev[n].xarr;
		my = dev[n].yarr;
		mz = dev[n].zarr;

		boundz<<< dim3(nx,ny,2), dim3(zpad,1,1), 0, dev[n].stream >>>(&dev[n].xa[dev[n].xbgn], dev[n].ya, dev[n].za, nz, mx, my, dev[n].C);
		sweepz<<< dim3(nz/z_zdiv,nx/z_xdiv,ny/z_ydiv), dim3(z_zthd,z_xdiv,z_ydiv), 2*sizeof(double)*z_zthd*z_xdiv*z_ydiv, dev[n].stream >>>
		      (dev[n].za, dev[n].zv, &dev[n].xa[dev[n].xbgn], &dev[n].xv[dev[n].xbgn], dev[n].ya, dev[n].yv, mx, my, dev[n].C, dev[n].F, dev[n].orb_shf, dt);
		update<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 0, dev[n].stream >>>
		      (dt, &dev[n].xa[dev[n].xbgn], dev[n].ya, dev[n].za, &dev[n].xv[dev[n].xbgn], dev[n].yv, dev[n].zv, mx, my, mz, dev[n].C, dev[n].F);
	}
	#endif

	#if ndim>1
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		mx = dev[n].xarr;
		my = dev[n].yarr;
		mz = dev[n].zarr;

		boundy<<< dim3(nx,nz,2), dim3(ypad,1,1), 0, dev[n].stream >>>(&dev[n].xa[dev[n].xbgn], dev[n].ya, dev[n].za, ny, mx, my, dev[n].C);
		sweepy<<< dim3(ny/y_ydiv,nx/y_xdiv,nz/y_zdiv), dim3(y_ythd,y_xdiv,y_zdiv), 2*sizeof(double)*y_ythd*y_xdiv*y_zdiv, dev[n].stream >>>
		      (dev[n].ya, dev[n].yv, &dev[n].xa[dev[n].xbgn], &dev[n].xv[dev[n].xbgn], dev[n].za, dev[n].zv, mx, my, dev[n].C, dev[n].F, dev[n].orb_rot, dt);

		update<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 0, dev[n].stream >>>
		      (dt, &dev[n].xa[dev[n].xbgn], dev[n].ya, dev[n].za, &dev[n].xv[dev[n].xbgn], dev[n].yv, dev[n].zv, mx, my, mz, dev[n].C, dev[n].F);
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

		mx = dev[n].xarr;
		my = dev[n].yarr;
		mz = dev[n].zarr;

		sweepx<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 2*sizeof(double)*x_xthd*x_ydiv*x_zdiv, dev[n].stream >>>
		      (dev[n], dt);

		update<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 0, dev[n].stream >>>
		      (dt, &dev[n].xa[dev[n].xbgn], dev[n].ya, dev[n].za, &dev[n].xv[dev[n].xbgn], dev[n].yv, dev[n].zv, mx, my, mz, dev[n].C, dev[n].F);
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

		mx = dev[n].xarr;
		my = dev[n].yarr;
		mz = dev[n].zarr;

		boundy<<< dim3(nx,nz,2), dim3(ypad,1,1), 0, dev[n].stream >>>(&dev[n].xa[dev[n].xbgn], dev[n].ya, dev[n].za, ny, mx, my, dev[n].C);
		advecy<<< dim3(ny/y_ydiv,nx/y_xdiv,nz/y_zdiv), dim3(y_ythd,y_xdiv,y_zdiv), 2*sizeof(double)*y_ythd*y_xdiv*y_zdiv, dev[n].stream >>>
		      (dev[n].ya, dev[n].yv, &dev[n].xa[dev[n].xbgn], &dev[n].xv[dev[n].xbgn], dev[n].za, dev[n].zv, mx, my, dev[n].C, dev[n].F, dev[n].orb_res, dt);

		updadv<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 0, dev[n].stream >>>
		      (dt, &dev[n].xa[dev[n].xbgn], dev[n].ya, dev[n].za, &dev[n].xv[dev[n].xbgn], dev[n].yv, dev[n].zv, mx, my, mz, dev[n].C, dev[n].F);
	}
	#endif
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

		planet_evo<<< 1, n_planet, 0, dev[n].stream >>> (dev[n].planets, time+hdt, dt);

		source<<< dim3((nx+1023)/1024,ny,nz), dim3(1024,1,1), 0, dev[n].stream >>>
		      (dev[n], hdt);

		#ifdef kill_flag
		killwave<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 0, dev[n].stream >>>
		        (&dev[n].xa[dev[n].xbgn], dev[n].ya, dev[n].za, mx, my, mz, dev[n].C, dt);
		#endif
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);

	return;
}

