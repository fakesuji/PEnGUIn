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
	int imax = blockDim.x;
	int i = threadIdx.x;

	flux.r *= gfac;
	flux.p *= gfac;
	flux.u *= gfac;
	flux.v *= gfac;
	flux.w *= gfac;

	//double old_flux = flux.r;
	double vol;

	if (i>=npad && i<imax-npad+1)
	{
		if (flux.r>0.0)
		{
			vol = 0.95*get_volume_dev(geom, xa[i-1], dx[i-1]);
			if (flux.r*dt>r[i-1]*vol)
			{
				//flux.r = r[i-1]*vol/dt;
				//flux.p *= flux.r/old_flux;
				//flux.u *= flux.r/old_flux;
				//flux.v *= flux.r/old_flux;
				//flux.w *= flux.r/old_flux;
				flux.r = 0.0;
				flux.p = 0.0;
				flux.u = 0.0;
				flux.v = 0.0;
				flux.w = 0.0;
				printf("Error: flux too high; consider lowering the Courant number.\n");
			}
		}
		else if (flux.r<0.0) 
		{
			vol = 0.95*get_volume_dev(geom, xa[i], dx[i]);
			if (-flux.r*dt>r[i]*vol)
			{
				//flux.r = -r[i]*vol/dt;
				//flux.p *= flux.r/old_flux;
				//flux.u *= flux.r/old_flux;
				//flux.v *= flux.r/old_flux;
				//flux.w *= flux.r/old_flux;
				flux.r = 0.0;
				flux.p = 0.0;
				flux.u = 0.0;
				flux.v = 0.0;
				flux.w = 0.0;
				printf("Error: flux too high; consider lowering the Courant number.\n");
			}
		}
	}	
	return;
}

__device__ void net_source(Cell &flux, int geom, double* r, double* xa, double* dx, double pres, double uprs, double dt)
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

	check_flux(flux, geom, r, xa, dx, gfac, dt);
	__syncthreads();

	flux.r = net_flux(tmp1,flux.r);
	__syncthreads();
	flux.p = net_flux(tmp1,flux.p) + de;
	__syncthreads();
	flux.u = net_flux(tmp1,flux.u) + du;
	__syncthreads();
	flux.v = net_flux(tmp1,flux.v);
	__syncthreads();
	flux.w = net_flux(tmp1,flux.w);
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

	#if recon_flag==2
	flatten(r, p, u);
	__syncthreads();
	#endif

	//p[i] = p[i]/r[i];
	//__syncthreads();

	if (i>=npad && i<imax+1-npad)
	{
		us = 0.5*dt*force;
		if (geom>2) dt /= rad;

		set_state(i, geom, xa, dx, dv, rad, r, p, u, v, w, dt, us, S);

		wave_speeds(S, pm, sl, sm, sr, us);

		//set_L_state_passive(i-1, geom, xa, dx, dv, rad, sl, sm, u, v, w, dt, S);
		//set_R_state_passive(  i, geom, xa, dx, dv, rad, sr, sm, u, v, w, dt, S);
//if (xa[i]<2.0 && xa[i]>1.993) printf("%f:\n %.10e, %.10e, %.10e, %.10e\n %.10e, %.10e, %.10e, %.10e\n %.10e, %.10e, %.10e\n", rad, S.rl, S.pl, S.ul, S.vl, S.rr, S.pr, S.ur, S.vr, sl, sm, sr);

		HLLC_fluxes(S, pm, sl, sm, sr, flux, pres, uprs);

		if (isnan(flux.r) && !isnan(u[i-1]*u[i]*u[i+1]*u[i-2]*p[i-1]*p[i]*p[i+1]*p[i-2]*r[i-1]*r[i]*r[i+1]*r[i-2]))
		{		
			printf("Error: flux nan, %e| %e, %e, %e, %e, %e|\n %e, %e, %e, %e|\n %e, %e, %e, %e|\n",
			pm,sl,S.ul,sm,S.ur,sr,
			r[i-2],r[i-1],r[i],r[i+1],
			p[i-2],p[i-1],p[i],p[i+1]);
			set_state(i, geom, xa, dx, dv, rad, r, p, u, v, w, dt, us, S, true);
		}
	}
	__syncthreads();

	net_source(flux, geom, r, xa, dx, pres, uprs, dt);

	return flux;
}

