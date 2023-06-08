#include <unistd.h>
#include <iostream>

#include "parameters.h"
#include "structs.h"
#include "init.h"
#include "boundary.h"

__device__ void get_deri_1D(double* u, double* r, double &dudx, double &d2udx2)
{
	dudx = (u[2]-u[0])/(r[2]-r[0]);
	d2udx2 = ((u[2]-u[1])/(r[2]-r[1]) - (u[1]-u[0])/(r[1]-r[0])) / (0.5*(r[2]-r[0]));
	return;
}

__device__ void get_deri_2D(double* u, double* r, double* t, double &dudy, double &d2udy2, double &d2udxy)
{
	dudy = (u[7]-u[4])/(t[2]-t[0])/r[1];
	d2udy2 = ((u[7]-u[1])/(t[2]-t[1]) - (u[1]-u[4])/(t[1]-t[0])) / (0.5*(t[2]-t[0])) / r[1] / r[1];
	d2udxy = (u[8] - u[6] - u[5] + u[3]) / (r[2]-r[0]) / (t[2]-t[0]) / r[1];
	return;
}

__device__ void get_deri_3D(double* u, double* r, double* t, double* z, double &dudz, double &d2udz2, double &d2udxz, double &d2udyz)
{
	dudz = (u[15]-u[10])/(z[2]-z[0]);
	d2udz2 = ((u[15]-u[1])/(z[2]-z[1]) - (u[1]-u[10])/(z[1]-z[0])) / (0.5*(z[2]-z[0]));
	d2udxz = (u[16] - u[14] - u[11] + u[9] ) / (z[2]-z[0]) / (r[2]-r[0]);
	d2udyz = (u[18] - u[17] - u[13] + u[12]) / (z[2]-z[0]) / (t[2]-t[0]) / r[1];
	return;
}


__global__ void viscous_forces(Grid G, Cell* C)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + 1;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;
	#if ndim>1
	j++;
	#endif
	#if ndim>2
	k++;
	#endif
	if (i>=G.xarr-1) return;

	double nu0, nu1, nu2;
	int s;
	#if ndim==1
	s = 3;
	double u[3], nr[3];
	double r[3];
	double dudx, d2udx2;
	double dndx;
	#elif ndim==2
	s = 9;
	double u[9], v[9], nr[9];
	double r[3], t[3];
	double dudx, d2udx2, dudy, d2udy2, d2udxy;
	double dvdx, d2vdx2, dvdy, d2vdy2, d2vdxy;
	double dndx, dndy;
	#elif ndim==3
	s = 19;
	double u[19], v[19], w[19], nr[19];
	double r[3], t[3], z[3];
	double dudx, d2udx2, dudy, d2udy2, d2udxy, dudz, d2udz2, d2udxz, d2udyz;
	double dvdx, d2vdx2, dvdy, d2vdy2, d2vdxy, dvdz, d2vdz2, d2vdxz, d2vdyz;
	double dwdx, d2wdx2, dwdy, d2wdy2, d2wdxy, dwdz, d2wdz2, d2wdxz, d2wdyz;
	double dndx, dndy, dndz;
	#endif

	int ind;
	for (int n=0; n<3; n++)
	{
		r[n] = G.get_xc(i+n-1);
		#if ndim>1
		t[n] = G.get_yc(j+n-1);
		#endif
		#if ndim>2
		z[n] = G.get_zc(k+n-1);
		#endif
	}

	nu0 = get_nu(r[0],G.get_yc(j),G.get_zc(k));
	nu1 = get_nu(r[1],G.get_yc(j),G.get_zc(k));
	nu2 = get_nu(r[2],G.get_yc(j),G.get_zc(k));

	#if ndim>2
	for (int n=0; n<3; n++)
	{
		z[n] *= r[1];
	}
	#endif

	for (int n=0; n<s; n++)
	{
		if      (n<3)  ind = G.get_ind(i+n-1,j,k);
		else if (n<6)  ind = G.get_ind(i+n-4,j-1,k);
		else if (n<9)  ind = G.get_ind(i+n-7,j+1,k);
		else if (n<12) ind = G.get_ind(i+n-10,j,k-1);
		else if (n<13) ind = G.get_ind(i,j-1,k-1);
		else if (n<14) ind = G.get_ind(i,j+1,k-1);
		else if (n<17) ind = G.get_ind(i+n-15,j,k+1);
		else if (n<18) ind = G.get_ind(i,j-1,k+1);
		else if (n<19) ind = G.get_ind(i,j+1,k+1);

		nr[n] = C[ind].r;
		u[n]  = C[ind].u;
		#if ndim>1
		v[n]  = C[ind].v;
		#endif
		#if ndim>2
		w[n]  = C[ind].w;
		#endif

		if (n==0 || n==3 || n==6 || n==9 || n==14)
		{
			nr[n] *= nu0;
		}
		else if (n==2 || n==5 || n==8 || n==11 || n==16)
		{
			nr[n] *= nu2;
		}
		else
		{
			nr[n] *= nu1;
		}	
	}

	get_deri_1D(u, r, dudx, d2udx2);
	dndx = (nr[2]-nr[0])/(r[2]-r[0]);
	#if ndim>1
	get_deri_1D(v, r, dvdx, d2vdx2);
	get_deri_2D(v, r, t, dvdy, d2vdy2, d2vdxy);
	get_deri_2D(u, r, t, dudy, d2udy2, d2udxy);
	dndy = (nr[7]-nr[4])/(t[2]-t[0])/r[1];
	#endif
	#if ndim>2
	get_deri_1D(w, r, dwdx, d2wdx2);
	get_deri_2D(w, r, t, dwdy, d2wdy2, d2wdxy);
	get_deri_3D(w, r, t, z, dwdz, d2wdz2, d2wdxz, d2wdyz);
	get_deri_3D(v, r, t, z, dvdz, d2vdz2, d2vdxz, d2vdyz);
	get_deri_3D(u, r, t, z, dudz, d2udz2, d2udxz, d2udyz);
	dndz = (nr[15]-nr[10])/(z[2]-z[0]);
	#endif

	////////////////////////////////////////////////////////////////////////////

	double sxx = 2.0*nr[1]*dudx;

	double fx;
	fx  = 2.0*dndx*dudx + 2.0*nr[1]*d2udx2 + sxx/r[1];

	////////////////////////////////////////////////////////////////////////////

	#if ndim>1
	double syy = 2.0*nr[1]*(dvdy+u[1]/r[1]);
	double sxy = nr[1]*(dudy + dvdx - v[1]/r[1]);

	double fy;
	fy  = 2.0*dndy*dvdy + 2.0*nr[1]*d2vdy2 + 2.0*nr[1]*dudy/r[1] + 2.0*dndy*u[1]/r[1];

	fx += dndy*(dudy + dvdx) + nr[1]*(d2udy2 + d2vdxy) - dndy*v[1]/r[1] - nr[1]*dvdy/r[1] - syy/r[1];
	fy += dndx*(dudy + dvdx) + nr[1]*(d2udxy + d2vdx2) - dndx*v[1]/r[1] - nr[1]*dvdx/r[1] + nr[1]*v[1]/r[1]/r[1] + 2.0*sxy/r[1];
	#endif

	////////////////////////////////////////////////////////////////////////////

	#if ndim>2
	double sxz = nr[1]*(dudz + dwdx);

	double fz;
	fz  = 0.0;//2.0*dndz*dwdz + 2.0*nr[1]*d2wdz2;

	//fx += dndz*(dudz + dwdx) + nr[1]*(d2udz2 + d2wdxz);
	//fz += dndx*(dudz + dwdx) + nr[1]*(d2udz2 + d2wdxz) + sxz/r[1];

	//fy += dndz*(dvdz + dwdy) + nr[1]*(d2vdz2 + d2wdyz); 
	//fz += dndy*(dvdz + dwdy) + nr[1]*(d2vdyz + d2wdy2);

	#endif

	////////////////////////////////////////////////////////////////////////////

	//if (i==100) printf("%i, %e: %e, %e, %e, %e, \n%e, %e, %e\n%e, %e, %e\n%e, %e, %e\n",i,r[1],dudx,dudy,dvdx,dvdy,  d2udxy,d2udx2,d2udy2,  d2vdxy,d2vdx2,d2vdy2, dndx,dndy,nr[1]);

	ind = G.get_ind(i,j,k);
	G.vis_tensor[ndim*ind+0] = fx;
	#if ndim>1
	G.vis_tensor[ndim*ind+1] = fy;
	#endif
	#if ndim>2
	G.vis_tensor[ndim*ind+2] = fz; 
	#endif

	return;
}

__global__ void viscous_tensor(Grid G, Cell* C)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + 1;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;
	#if ndim>1
	j++;
	#endif
	#if ndim>2
	k++;
	#endif
	if (i>=G.xarr-1) return;

	int s;
	#if ndim==1
	s = 1;
	#elif ndim==2
	s = 3;
	#elif ndim==3
	s = 6;
	#endif

	int ind;
	double up, vp, wp;
	double um, vm, wm;

	double r = G.get_xc(i);
	double dr = G.get_xc(i+1)-G.get_xc(i-1);

	ind = G.get_ind(i+1, j, k);
	up = C[ind].u;
	vp = C[ind].v;
	wp = C[ind].w;

	ind = G.get_ind(i-1, j, k);
	um = C[ind].u;
	vm = C[ind].v;
	wm = C[ind].w;

	double dudx = (up-um)/dr;
	#if ndim>1
	double dvdx = (vp-vm)/dr;
	#endif
	#if ndim>2
	double dwdx = (wp-wm)/dr;
	#endif

	#if ndim>1
	double dphi = G.get_yc(j+1)-G.get_yc(j-1);

	ind = G.get_ind(i, j+1, k);
	up = C[ind].u;
	vp = C[ind].v;
	wp = C[ind].w;

	ind = G.get_ind(i, j-1, k);
	um = C[ind].u;
	vm = C[ind].v;
	wm = C[ind].w;

	double dudy = (up-um)/dphi/r;
	double dvdy = (vp-vm)/dphi/r;
	#if ndim>2
	double dwdy = (wp-wm)/dphi/r;
	#endif
	#endif

	#if ndim>2
	double dz   = G.get_zc(k+1)-G.get_zc(k-1);

	ind = G.get_ind(i, j, k+1);
	up = C[ind].u;
	vp = C[ind].v;
	wp = C[ind].w;

	ind = G.get_ind(i, j, k-1);
	um = C[ind].u;
	vm = C[ind].v;
	wm = C[ind].w;

	double dudz = (up-um)/dz/r;
	double dvdz = (vp-vm)/dz/r;
	double dwdz = (wp-wm)/dz/r;
	#endif

	ind = G.get_ind(i, j, k);
	double rho = C[ind].r;
	double u = C[ind].u;
	double v = C[ind].v;

	double nr = get_nu(r,G.get_yc(j),G.get_zc(k))*rho;

	G.vis_tensor[s*ind+0] = nr*(dudx);

	#if ndim>1
	G.vis_tensor[s*ind+1] = nr*(dudy + dvdx - v/r);
	G.vis_tensor[s*ind+2] = nr*(dvdy + u/r);
	#endif

	#if ndim>2
	G.vis_tensor[s*ind+3] = nr*(dudz + dwdx);
	G.vis_tensor[s*ind+4] = nr*(dvdz + dwdy);
	G.vis_tensor[s*ind+5] = nr*(dwdz);
	#endif

	//printf("%f,%e,%e\n",r,nr,dudx);
	return;
}

/*
__device__ double viscous_fx(Grid G, Cell* C, int i, int j, int k)
{

	int s;
	#if ndim==1
	s = 1;
	#elif ndim==2
	s = 3;
	#elif ndim==3
	s = 6;
	#endif

	int ind;
	double fx;

	ind = G.get_ind(i+1, j, k);
	double* tf = &G.vis_tensor[s*ind];

	ind = G.get_ind(i-1, j, k);
	double* tb = &G.vis_tensor[s*ind];

	#if ndim>1
	ind = G.get_ind(i, j+1, k);
	double* tr = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j-1, k);
	double* tl = &G.vis_tensor[s*ind];
	#endif

	#if ndim>2
	ind = G.get_ind(i, j, k+1);
	double* tt = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j, k-1);
	double* td = &G.vis_tensor[s*ind];
	#endif

	ind = G.get_ind(i, j, k);
	double* t0 = &G.vis_tensor[s*ind];

	double r = G.get_xc(i);
	double dr   = G.get_xc(i+1)-G.get_xc(i-1);
	fx = 2.0*t0[0]/r + 2.0*(tf[0]-tb[0])/dr;

	#if ndim>1
	double dphi = G.get_yc(j+1)-G.get_yc(j-1);
	fx += (tr[1]-tl[1])/dphi/r - 2.0*t0[2]/r;
	#endif

	#if ndim>2
	double dz   = G.get_zc(k+1)-G.get_zc(k-1);
	fx += (tt[3]-td[3])/dz;
	#endif

	return fx;
}

__device__ double viscous_fy(Grid G, Cell* C, int i, int j, int k)
{
	int s;
	#if ndim==1
	s = 1;
	#elif ndim==2
	s = 3;
	#elif ndim==3
	s = 6;
	#endif

	int ind;
	double fy;

	ind = G.get_ind(i+1, j, k);
	double* tf = &G.vis_tensor[s*ind];

	ind = G.get_ind(i-1, j, k);
	double* tb = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j+1, k);
	double* tr = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j-1, k);
	double* tl = &G.vis_tensor[s*ind];

	#if ndim>2
	ind = G.get_ind(i, j, k+1);
	double* tt = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j, k-1);
	double* td = &G.vis_tensor[s*ind];
	#endif

	ind = G.get_ind(i, j, k);
	double* t0 = &G.vis_tensor[s*ind];

	double r = G.get_xc(i);
	double dr   = G.get_xc(i+1)-G.get_xc(i-1);
	double dphi = G.get_yc(j+1)-G.get_yc(j-1);

	fy = 2.0*t0[1]/r + (tf[1]-tb[1])/dr + 2.0*(tr[2]-tl[2])/dphi/r;

	#if ndim>2
	double dz   = G.get_zc(k+1)-G.get_zc(k-1);
	fy += (tt[4]-td[4])/dz;
	#endif

	return fy;
}

__device__ double viscous_fz(Grid G, Cell* C, int i, int j, int k)
{
	int s;
	#if ndim==1
	s = 1;
	#elif ndim==2
	s = 3;
	#elif ndim==3
	s = 6;
	#endif

	int ind;
	double fz;

	ind = G.get_ind(i+1, j, k);
	double* tf = &G.vis_tensor[s*ind];

	ind = G.get_ind(i-1, j, k);
	double* tb = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j+1, k);
	double* tr = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j-1, k);
	double* tl = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j, k+1);
	double* tt = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j, k-1);
	double* td = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j, k);
	double* t0 = &G.vis_tensor[s*ind];

	double r = G.get_xc(i);
	double dr   = G.get_xc(i+1)-G.get_xc(i-1);
	double dphi = G.get_yc(j+1)-G.get_yc(j-1);
	double dz   = G.get_zc(k+1)-G.get_zc(k-1);

	fz = t0[3]/r + (tf[3]-tb[3])/dr + (tr[4]-tl[4])/dphi/r + 2.0*(tt[5]-td[5])/dz;

	return fz;
}
*/

__device__ double viscous_heat(Grid G, Cell* C, int i, int j, int k)
{
	int s;
	#if ndim==1
	s = 1;
	#elif ndim==2
	s = 3;
	#elif ndim==3
	s = 6;
	#endif

	int ind;
	double H;

	ind = G.get_ind(i, j, k);
	double* T = &G.vis_tensor[s*ind];

	double rho = C[ind].r;
	double nr = get_nu(G.get_xc(i),G.get_yc(j),G.get_zc(k))*rho;

	H = 2.0*T[0]*T[0] + 2.0*T[2]*T[2] + 2.0*T[5]*T[5] + T[1]*T[1] + T[3]*T[3] + T[4]*T[4];
	H /= nr;

	return H;
}

void viscosity_tensor_evaluation1(Grid* dev)
{
	int nx, ny, nz;

	#if ndim>2
	boundz(dev);
	#endif
	#if ndim>1
	boundy(dev);
	#endif
	boundx(dev);
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xarr-2;
		#if ndim>1
		ny = dev[n].yarr-2;
		#else
		ny = 1;
		#endif
		#if ndim>2
		nz = dev[n].zarr-2;
		#else
		nz = 1;
		#endif
		viscous_forces<<< dim3((nx+255)/256,ny,nz), 256, 0, dev[n].stream >>>(dev[n], dev[n].C);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
}

void viscosity_tensor_evaluation2(Grid* dev)
{
	int nx, ny, nz;

	#if ndim>2
	boundz2(dev);
	#endif
	#if ndim>1
	boundy2(dev);
	#endif
	boundx2(dev);
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xarr-2;
		#if ndim>1
		ny = dev[n].yarr-2;
		#else
		ny = 1;
		#endif
		#if ndim>2
		nz = dev[n].zarr-2;
		#else
		nz = 1;
		#endif
		viscous_forces<<< dim3((nx+255)/256,ny,nz), 256, 0, dev[n].stream >>>(dev[n], dev[n].T);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
}

