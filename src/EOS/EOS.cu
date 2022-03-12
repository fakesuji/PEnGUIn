#include "parameters.h"
#include "structs.h"

__host__ __device__ double get_energy(double r, double p, double u, double v, double w)
{
	double e;
	#if EOS_flag==2
	e = p/(gamm*r);
	#if internal_e_flag==0
	e+= u*u/2.0 + v*v/2.0 + w*w/2.0;
	#endif
	#elif EOS_flag==1
	e = p/pow(r,gam);
	#else
	e = p/r;
	#endif
	return e;
}
