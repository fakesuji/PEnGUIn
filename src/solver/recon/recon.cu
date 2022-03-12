#include "PEM.cu"
#include "PLM.cu"

__device__ void get_CON_parameters(int i, int geom, double* x, double* dx, double* dv, double* a, double* par)
{
	#if recon_flag==0
	get_PEM_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==1
	get_PLM_parameters(i, geom, x, dx, dv, a, par);
	#endif
	return;
}

__device__ double get_CON_aveR(int geom, double rL, double r0, double rR, double* par)
{
	if (r0<rL) {printf("error! %f, %f, %f\n", rL, r0, rR); r0 = rL;}
	if (r0>rR) {printf("error! %f, %f, %f\n", rL, r0, rR); r0 = rR;}
	#if recon_flag==0
	return get_PEM_aveR(geom, rL, r0, rR, par);
	#elif recon_flag==1
	return get_PLM_aveR(geom, rL, r0, rR, par);
	#endif
}

__device__ double get_CON_aveL(int geom, double rL, double r0, double rR, double* par)
{
	if (r0<rL) {printf("error! %f, %f, %f\n", rL, r0, rR); r0 = rL;}
	if (r0>rR) {printf("error! %f, %f, %f\n", rL, r0, rR); r0 = rR;}
	#if recon_flag==0
	return get_PEM_aveL(geom, rL, r0, rR, par);
	#elif recon_flag==1
	return get_PLM_aveL(geom, rL, r0, rR, par);
	#endif
}

__device__ void get_PRM_parameters(int i, int geom, double* x, double* dx, double* dv, double* a, double* par)
{
	#if recon_flag==0
	get_PLM_parameters(i, geom, x, dx, dv, a, par);
	#elif recon_flag==1
	get_PLM_parameters(i, geom, x, dx, dv, a, par);
	#endif
}

__device__ double get_PRM_aveR(int geom, double rL, double r0, double rR, double* par)
{
	#if recon_flag==0
	return get_PLM_aveR(geom, rL, r0, rR, par);
	#elif recon_flag==1
	return get_PLM_aveR(geom, rL, r0, rR, par);
	#endif
}

__device__ double get_PRM_aveL(int geom, double rL, double r0, double rR, double* par)
{
	#if recon_flag==0
	return get_PLM_aveL(geom, rL, r0, rR, par);
	#elif recon_flag==1
	return get_PLM_aveL(geom, rL, r0, rR, par);
	#endif
}
