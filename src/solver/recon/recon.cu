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

__device__ double get_CON_aveR(int geom, double x, double* par, double lx, double ly)
{
	#if recon_flag==0
	return get_PEM_aveR(geom, x, par, lx, ly);
	#elif recon_flag==1
	return get_PLM_aveR(geom, x, par);
	#elif recon_flag==2
	return get_PPM_aveR(geom, x, par);
	#endif
}

__device__ double get_CON_aveL(int geom, double x, double* par, double lx, double ly)
{
	#if recon_flag==0
	return get_PEM_aveL(geom, x, par, lx, ly);
	#elif recon_flag==1
	return get_PLM_aveL(geom, x, par);
	#elif recon_flag==2
	return get_PPM_aveL(geom, x, par);
	#endif
}

__device__ void get_PRM_parameters(int i, int geom, double* x, double* dx, double* dv, double* a, double* par)
{
	#if recon_flag==0
	get_PEM_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==1
	get_MOC_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==2
	get_PPM_parameters(i, geom, x, dx, dv, a, par);
	#endif
}

__device__ double get_PRM_aveR(int geom, double x, double* par, double lx, double ly)
{
	#if recon_flag==0
	return get_PEM_aveR(geom, x, par, lx, ly);
	#elif recon_flag==1
	return get_PLM_aveR(geom, x, par);
	#elif recon_flag==2
	return get_PPM_aveR(geom, x, par);
	#endif
}

__device__ double get_PRM_aveL(int geom, double x, double* par, double lx, double ly)
{
	#if recon_flag==0
	return get_PEM_aveL(geom, x, par, lx, ly);
	#elif recon_flag==1
	return get_PLM_aveL(geom, x, par);
	#elif recon_flag==2
	return get_PPM_aveL(geom, x, par);
	#endif
}

__device__ double get_slopeR(int geom, double x1, double x2, double* par)
{
	#if recon_flag==0
	return get_PEM_slopeR(geom, x1, x2, par);
	#elif recon_flag==1
	return get_PLM_slopeR(geom, x1, x2, par);
	#elif recon_flag==2
	return get_PPM_slopeR(geom, x1, x2, par);
	#endif
}

__device__ double get_slopeL(int geom, double x1, double x2, double* par)
{
	#if recon_flag==0
	return get_PEM_slopeL(geom, x1, x2, par);
	#elif recon_flag==1
	return get_PLM_slopeL(geom, x1, x2, par);
	#elif recon_flag==2
	return get_PPM_slopeL(geom, x1, x2, par);
	#endif
}
