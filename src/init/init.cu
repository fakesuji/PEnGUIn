#include "unistd.h"
#include "iostream"
#include "parameters.h"
#include "structs.h"
#include "geom.h"

__host__ __device__ double get_cs2(double x, double y, double z)
{
	double coeff = sc_h*sc_h*(1.0/planet_radius);
	#if geomx == 0
	return coeff;
	#elif geomx == 1
	return coeff*pow(x/planet_radius,-p_beta);
	#elif geomx == 2
	return coeff*pow(x*sin(z)/planet_radius,-p_beta);
	#endif
}


__host__ __device__ double get_h(double x, double y, double z)
{
	#if geomx == 0
	return sc_h;
	#elif geomx == 1
	return sqrt(get_cs2(x, y, z)*pow(x,3));
	#elif geomx == 2
	return sqrt(get_cs2(x, y, z)*pow(x*sin(z),3));
	#endif
}

__host__ __device__ double get_nu(double x, double y, double z)
{
	#ifdef visc_flag

	double coeff = ss_alpha*sqrt(get_cs2(planet_radius, y, z))*get_h(planet_radius, y, z);
	#if geomx == 0
	return coeff;
	#elif geomx == 1
	return coeff*pow(x/planet_radius,p_alpha);//pow(x,-p_beta+1.5);
	#elif geomx == 2
	double rad_cyl = x*sin(z);
	return coeff*pow(rad_cyl/planet_radius,p_alpha);//pow(rad_cyl,-p_beta+1.5);
	#endif

	#else

	return 0.0;

	#endif
}

__host__ __device__ double get_vertical(double x, double y, double z)
{
	#if ndim>2
		#if geomx == 0
		return exp(-z*z/(2.0*sc_h*sc_h));
		#elif geomx == 1
		double cs2 = get_cs2(x,y,z);
		double s = sqrt(x*x + z*z);
		return exp((1.0/s - 1.0/x)/cs2);
		#elif geomx == 2
		double cs2 = get_cs2(x,y,z);
		double r = x * sin(z);
		return exp((1.0/x - 1.0/r)/cs2);
		#endif
	#else
		return 1.0;
	#endif

}

__host__ __device__ double get_r(double x, double y, double z)
{
	#if init_flag == 0
	if (x>=0.3 && x <=0.7) return 2.0;
	else                   return 1.0;

	#elif init_flag == 1
	return 1.0;

	#elif init_flag == 2
	return pow(x/planet_radius,-p_alpha)*get_vertical(x,y,z);

	#elif init_flag == 3
	return pow(x*sin(z)/planet_radius,-p_alpha)*get_vertical(x,y,z);

	#elif init_flag == 4
	return 1.0;

	#elif init_flag == 5

	#elif init_flag == 6

	#elif init_flag == 7

	#elif init_flag == 8

	#endif
}


__host__ __device__ double get_p(double x, double y, double z)
{
	#if init_flag == 0
	return 1.0;

	#elif init_flag == 1
	#ifdef rev_flag
	x = 1.0-x;
	#endif
	if      (x<=0.1) return 100.0;
	else if (x>=0.9) return 1000.0;
	else             return 0.01;

	#elif init_flag == 2
	return get_cs2(x,y,z)*get_r(x,y,z);

	#elif init_flag == 3
	return get_cs2(x,y,z)*get_r(x,y,z);

	#elif init_flag == 4
	if (x<0.5 && y<0.5) return 2.0;
	else return 1.0;

	#elif init_flag == 5

	#elif init_flag == 6

	#elif init_flag == 7

	#elif init_flag == 8

	#endif
}

__host__ __device__ double get_u(double x, double y, double z)
{
	#if init_flag == 0
	return 1.0;

	#elif init_flag == 1
	return 0.0;

	#elif init_flag == 2
	return -1.5*get_nu(x,y,z)/x;

	#elif init_flag == 3
	double r = x * sin(z);
	return -1.5*get_nu(x,y,z)/r;

	#elif init_flag == 4
	return 0.0;

	#elif init_flag == 5

	#elif init_flag == 6

	#elif init_flag == 7

	#elif init_flag == 8

	#endif
}

__host__ __device__ double get_v(double x, double y, double z)
{
	#if init_flag == 0
	return 0.0;

	#elif init_flag == 1
	return 0.0;

	#elif init_flag == 2
	double rho = get_r(x,y,z);
	double dr = 1.1514178e-7;
	double dP_dr = (get_p(x+dr,y,z)-get_p(x-dr,y,z))/(2.0*dr);
	return sqrt(1.0/x + (x/rho)*dP_dr);

	#elif init_flag == 3
	double rho = get_r(x,y,z);
	double dr = 1.1514178e-7;
	double dP_dr = (get_p(x+dr,y,z)-get_p(x-dr,y,z))/(2.0*dr);
	return sqrt(1.0/x + (x/rho)*dP_dr);

	#elif init_flag == 4
	return 0.0;

	#elif init_flag == 5

	#elif init_flag == 6

	#elif init_flag == 7

	#elif init_flag == 8

	#endif
}


__host__ __device__ double get_w(double x, double y, double z)
{
	#if init_flag == 0
	return 0.0;

	#elif init_flag == 1
	return 0.0;

	#elif init_flag == 2
	return 0.0;

	#elif init_flag == 3
	return 0.0;

	#elif init_flag == 4
	return 0.0;

	#elif init_flag == 5

	#elif init_flag == 6

	#elif init_flag == 7

	#elif init_flag == 8

	#endif
}

__host__ __device__ Cell init_C(double x, double y, double z)
{
	Cell Q;

	Q.r = get_r(x,y,z);
	Q.p = get_p(x,y,z);
	Q.u = get_u(x,y,z);
	Q.v = get_v(x,y,z);
	Q.w = get_w(x,y,z);

	return Q;
}

__host__ __device__ Cell init_C(double x0, double x1, double y0, double y1, double z0, double z1)
{
	return init_C(0.5*(x0+x1),0.5*(y0+y1),0.5*(z0+z1));
}

void make_grid(double* a, double* v, double amin, double amax, int res, int pad, int geom, int grid)
{
	if      (grid==0) linspace(&a[pad], amin, amax, res+1);
	else if (grid==1) logspace(&a[pad], amin, amax, res+1); 

	for (int i = 0; i<pad; i++)
	{
		if (grid==1)
		{
			a[i] = amin*exp(log(amax/amin)*(double)(i-pad)/(double)(res));
			a[res+pad+1+i] = amin*exp(log(amax/amin)*(double)(res+1+i)/(double)(res));
		}
		else
		{
			a[i] = a[pad] - (a[2*pad-i]-a[pad]);
			a[res+pad+1+i] = a[res+pad] + (a[res+pad]-a[res+pad-1-i]);
		}
	}

	for (int i=0; i<res+2*pad; i++)
	{
		v[i] = get_volume(geom,a[i],a[i+1]-a[i]);
		//printf("%i, %f,%f,%f\n",i,a[i],a[i+1],v[i]);
	}
	return;
}

void fill_grid(double* xa, double* xv, double* ya, double* yv, double* za, double* zv, int mx, int my, int mz, Cell* C)
{
	int ind;
	for (int i=0; i<mx; i++)
	for (int j=0; j<my; j++)
	for (int k=0; k<mz; k++)
	{
		ind = i + mx*(j + my*k);
		C[ind] = init_C(xa[i], xa[i+1], ya[j], ya[j+1], za[k], za[k+1]);
	}
	return;
}

void init(Grid* dev)
{
	for (int n=0; n<ndev; n++)
	{
		make_grid(dev[n].xa, dev[n].xv, xmin, xmax, xres, xpad, geomx, gridx);
		make_grid(dev[n].ya, dev[n].yv, ymin, ymax, yres, ypad, geomy, gridy);
		make_grid(dev[n].za, dev[n].zv, zmin, zmax, zres, zpad, geomz, gridz);
		fill_grid(&dev[n].xa[dev[n].xbgn], &dev[n].xv[dev[n].xbgn], 
		          &dev[n].ya[dev[n].ybgn], &dev[n].yv[dev[n].ybgn], 
		          &dev[n].za[dev[n].zbgn], &dev[n].zv[dev[n].zbgn], 
		          dev[n].xarr, dev[n].yarr, dev[n].zarr, dev[n].C);
	}

	return;
}
