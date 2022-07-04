#include "PEM.cu"
#include "PLM.cu"
#include "PPM.cu"

__device__ void get_CON_parameters(int i, int geom, double* x, double* dx, double* dv, double* a, double* par)
{
	#if recon_flag==0
	get_PEM_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==1
	get_MOC_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==2
	get_PPM_parameters(i, geom, x, dx, dv, a, par);
	#endif
	return;
}

__device__ double get_CON_aveR(int geom, double rL, double r0, double rR, double* par)
{
	#if recon_flag==0
	return get_PEM_aveR(geom, rL, r0, rR, par);
	#elif recon_flag==1
	return get_PLM_aveR(geom, rL, r0, rR, par);
	#elif recon_flag==2
	return get_PPM_aveR(geom, rL, r0, rR, par);
	#endif
}

__device__ double get_CON_aveL(int geom, double rL, double r0, double rR, double* par)
{
	#if recon_flag==0
	return get_PEM_aveL(geom, rL, r0, rR, par);
	#elif recon_flag==1
	return get_PLM_aveL(geom, rL, r0, rR, par);
	#elif recon_flag==2
	return get_PPM_aveL(geom, rL, r0, rR, par);
	#endif
}

__device__ void get_PRM_parameters(int i, int geom, double* x, double* dx, double* dv, double* a, double* par)
{
	#if recon_flag==0
	get_VAN_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==1
	get_VAN_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==2
	get_PPM_parameters(i, geom, x, dx, dv, a, par);
	#endif
}

__device__ double get_PRM_aveR(int geom, double rL, double r0, double rR, double* par)
{
	#if recon_flag==0
	return get_PLM_aveR(geom, rL, r0, rR, par);
	#elif recon_flag==1
	return get_PLM_aveR(geom, rL, r0, rR, par);
	#elif recon_flag==2
	return get_PPM_aveR(geom, rL, r0, rR, par);
	#endif
}

__device__ double get_PRM_aveL(int geom, double rL, double r0, double rR, double* par)
{
	#if recon_flag==0
	return get_PLM_aveL(geom, rL, r0, rR, par);
	#elif recon_flag==1
	return get_PLM_aveL(geom, rL, r0, rR, par);
	#elif recon_flag==2
	return get_PPM_aveL(geom, rL, r0, rR, par);
	#endif
}

__device__ double get_distance_traveled_L(double* u_par, double dt)
{
	double m = u_par[2]-u_par[0];
	double v = u_par[0];

	if (v>0.0) return 0.0;
	else       
	{
		if (fabs(m)<0.0001)
		{
			return v*dt*( 1.0 - (m*dt/2.0) + (m*m*dt*dt/6.0) - (m*m*m*dt*dt*dt/24.0) );
		}
		else
		{
			return (v/m)*(exp(-m*dt)-1.0);
		}
	}
}

__device__ double get_distance_traveled_R(double* u_par, double dt)
{
	double m = u_par[0]-u_par[2];
	double v = u_par[2];

	if (v>0.0) 
	{
		if (fabs(m)<0.0001)
		{
			return v*dt*( 1.0 + (m*dt/2.0) + (m*m*dt*dt/6.0) + (m*m*m*dt*dt*dt/24.0) );
		}
		else
		{
			return (v/m)*(exp(m*dt)-1.0);
		}
	}
	else       return 0.0;
}
