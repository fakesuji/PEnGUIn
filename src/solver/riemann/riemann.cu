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

__device__ void net_source(Cell &flux, int geom, double* xa, double pres, double uprs)
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
                        double* r, double* p, double* u, double* v, double* w, double dt)
{
	int imax = blockDim.x;
	int i = threadIdx.x;

	State S;
	double pm, sl, sm, sr;
	Cell flux;
	double pres, uprs;

	if (i>=npad && i<imax+1-npad)
	{
		if (geom>2) dt /= rad;

		set_L_state(i-1, geom, xa, dx, dv, rad, r, p, u, v, w, dt, S);
		set_R_state(  i, geom, xa, dx, dv, rad, r, p, u, v, w, dt, S);
		wave_speeds(S, pm, sl, sm, sr);
		set_L_state_passive(i-1, geom, xa, dx, dv, rad, sl, sm, u, v, w, dt, S);
		set_R_state_passive(  i, geom, xa, dx, dv, rad, sr, sm, u, v, w, dt, S);

//if (xa[i]<2.0 && xa[i]>1.993) printf("%f:\n %.10e, %.10e, %.10e, %.10e\n %.10e, %.10e, %.10e, %.10e\n %.10e, %.10e, %.10e\n", rad, S.rl, S.pl, S.ul, S.vl, S.rr, S.pr, S.ur, S.vr, sl, sm, sr);

		HLLC_fluxes(S, pm, sl, sm, sr, flux, pres, uprs);
	}

	net_source(flux, geom, xa, pres, uprs);

	return flux;
}

