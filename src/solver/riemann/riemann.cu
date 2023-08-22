#include "states.cu"
#include "wave_speeds.cu"
#include "HLLC.cu"

__device__ double eval_gfac(int geom, double x)
{
	if      (geom==1) return x;
	else if (geom==2) return x*x;
	else if (geom==5) return sin(x);
	else              return 1.0;
}

//=======================================================================================

__device__ double net_flux(double* fluxes, double flux)
{
	int imax = blockDim.x;
	int i = threadIdx.x;
	int is = i + imax*threadIdx.y;

	fluxes[is] = flux;
	__syncthreads();
	if (i>=npad && i<imax-1) 
	{
		return fluxes[is]-fluxes[is+1];
	}
	else return 0.0;
}

__device__ void check_flux(Cell &flux, int geom, double* r, double* xa, double* dx, double gfac, double dt)
{
	flux.r *= gfac;
	flux.p *= gfac;
	flux.u *= gfac;
	flux.v *= gfac;
	flux.w *= gfac;
	return;
}

__device__ void net_source(Cell &flux, int geom, double* r, double* xa, double* dx, double pres, double uprs)
{
	int imax = blockDim.x;
	int i = threadIdx.x;
	int is = i + imax*threadIdx.y;

	double du=0.0, de=0.0;
	extern __shared__ double share[];
	double* tmp1 = &share[0];
	double* tmp2 = &share[blockDim.x*blockDim.y*blockDim.z];
	double gfac = eval_gfac(geom,xa[i]);

	if (i>=npad && i<imax+1-npad)
	{
		tmp1[is] = pres;
		tmp2[is] = gfac;
	}
	__syncthreads();

	if (i>=npad && i<imax-npad)
	{
		du = 0.5*(tmp2[is]+tmp2[is+1])*(tmp1[is]-tmp1[is+1]);
	}
	__syncthreads();

	#if EOS_flag==2
	if (i>=npad && i<imax+1-npad) tmp2[is] = uprs*gfac;
	__syncthreads();
	#if internal_e_flag==0
	if (i>=npad && i<imax-npad) de = tmp2[is]-tmp2[is+1];
	#elif internal_e_flag==1
	if (i>=npad && i<imax-npad) de = 0.5*(tmp1[is]+tmp1[is+1])*(tmp2[is]-tmp2[is+1]);
  	#endif
	#endif
	__syncthreads();

	flux.r = net_flux(tmp1,flux.r*gfac);
	__syncthreads();
	flux.p = net_flux(tmp1,flux.p*gfac) + de;
	__syncthreads();
	flux.u = net_flux(tmp1,flux.u*gfac) + du;
	__syncthreads();
	flux.v = net_flux(tmp1,flux.v*gfac);
	__syncthreads();
	flux.w = net_flux(tmp1,flux.w*gfac);
	__syncthreads();
	return;
}

//=======================================================================================

__device__ Cell riemann(int geom, double* xa, double* dx, double* dv, double rad,
                        double* r, double* p, double* u, double* v, double* w, double force, double dt)
{
	int imax = blockDim.x;
	int i = threadIdx.x;

	State S;
	double pm, sl, sm, sr, us;
	Cell flux;
	double pres, uprs;

	#if recon_flag==6
	flatten(r, p, u);
	__syncthreads();
	#endif

	if (i>=npad && i<imax+1-npad)
	{
		us = 0.5*dt*force;

		set_state(i, geom, xa, dx, dv, rad, r, p, u, v, w, dt/rad, us, S);
		wave_speeds(S, pm, sl, sm, sr, us);
		HLLC_fluxes(S, pm, sl, sm, sr, flux, pres, uprs);

		if (isnan(flux.r) && !isnan(u[i-1]*u[i]*u[i+1]*u[i-2]*p[i-1]*p[i]*p[i+1]*p[i-2]*r[i-1]*r[i]*r[i+1]*r[i-2]))
		{
			#ifndef silence_flag		
			printf("Error: flux nan, %e| %e, %e, %e, %e, %e|\n %e, %e, %e, %e|\n %e, %e, %e, %e|\n",
			pm,sl,S.ul,sm,S.ur,sr,
			r[i-2],r[i-1],r[i],r[i+1],
			p[i-2],p[i-1],p[i],p[i+1]);
			#endif
			set_state(i, geom, xa, dx, dv, rad, r, p, u, v, w, dt/rad, us, S, true);
		}
	}
	__syncthreads();

	net_source(flux, geom, r, xa, dx, pres, uprs);

	return flux;
}

