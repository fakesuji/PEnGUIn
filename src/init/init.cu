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
	if (x>=0.25 && x<=0.75) return 2.0;
	return 1.0;

	#elif init_flag == 6
	double A = 1.0e-6*sin(twopi*x);
	return 1.0 + A;

	#elif init_flag == 7
	double r1 = 1.0;
	double r2 = 2.0;
	double w = 0.002;

	if (y<=0.5)
	{
		return r1*(erf((0.25-y)/w)+1.0)/2.0 + r2*(erf((y-0.25)/w)+1.0)/2.0;
	}
	else
	{
		return r1*(erf((y-0.75)/w)+1.0)/2.0 + r2*(erf((0.75-y)/w)+1.0)/2.0;
	}

	#elif init_flag == 8
	return 1.0;

	#elif init_flag == 9
	if (x<0.125) return 3.857143;
	else       return 1.0 + 0.2*sin(twopi*x*8.0);

	#elif init_flag == 10
	if (x<1.0/6.0 + 0.57735026919*y) return 8.0;
	else                             return 1.4;

	#elif init_flag == 11
	if (x+y<0.15) return 0.125;
	else          return 1.0;

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
	if      (x<=0.1) return 1000.0;
	else if (x>=0.9) return 100.0;
	else             return 0.01;

	#elif init_flag == 2
	return get_cs2(x,y,z)*get_r(x,y,z);

	#elif init_flag == 3
	return get_cs2(x,y,z)*get_r(x,y,z);

	#elif init_flag == 4
	if (x<0.5 && y<0.5) return 2.0;
	else return 1.0;

	#elif init_flag == 5
	return 1.0;

	#elif init_flag == 6
	double A = 1.0e-6*sin(twopi*x);
	return (1.0 + gam*A)*(1.0 + A)/gam;

	#elif init_flag == 7
	return 2.5;

	#elif init_flag == 8
	return 3.0*gamm;

	#elif init_flag == 9
	if (x<0.125) return 10.33333;
	else       return 1.0;

	#elif init_flag == 10
	if (x<1.0/6.0 + 0.57735026919*y) return 116.5;
	else                             return 1.0;

	#elif init_flag == 11
	if (x+y<0.15) return 0.14;
	else          return 1.0;

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
	return 1.0;

	#elif init_flag == 6
	double A = 1.0e-6*sin(twopi*x);
	return A;

	#elif init_flag == 7
	double u1 = 0.5;
	double u2 =-0.5;
	double w = 0.002;

	if (y<=0.5)
	{
		return u1*(erf((0.25-y)/w)+1.0)/2.0 + u2*(erf((y-0.25)/w)+1.0)/2.0;
	}
	else
	{
		return u1*(erf((y-0.75)/w)+1.0)/2.0 + u2*(erf((0.75-y)/w)+1.0)/2.0;
	}


	#elif init_flag == 8
	if (x<=0.5) return -2.0;
	else        return 2.0;

	#elif init_flag == 9
	if (x<0.125) return 2.629369;
	else       return 0.0;

	#elif init_flag == 10
	if (x<1.0/6.0 + 0.57735026919*y) return 8.25*0.86602540378;
	else                             return 0.0;

	#elif init_flag == 11
	return 0.0;

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
	return 0.0;

	#elif init_flag == 6
	return 0.0;

	#elif init_flag == 7
/*
	double A = 1.0e-12;
	double v = 0.0;
	double k;

	for (int n=1; n<31; n++)
	{
		k = twopi*(double)n;
		if (y<0.25)
			//return A*((2.0*rand()/RAND_MAX)-1.0)*(exp(-k*(y+0.25))+exp(-k*(0.25-y)));
			v += k*sin(k*x)*(exp(-k*(y+0.25))+exp(-k*(0.25-y)));
		else if (y<0.75)
			//return A*((2.0*rand()/RAND_MAX)-1.0)*(exp(-k*(y-0.25))+exp(-k*(0.75-y)));
			v += k*sin(k*x)*(exp(-k*(y-0.25))+exp(-k*(0.75-y)));
		else
			//return A*((2.0*rand()/RAND_MAX)-1.0)*(exp(-k*(y-0.75))+exp(-k*(1.25-y)));
			v += k*sin(k*x)*(exp(-k*(y-0.75))+exp(-k*(1.25-y)));
	}
	return A*v;
*/
	return 1.0e-10*sin(twopi*x);
	#elif init_flag == 8
	return 0.0;

	#elif init_flag == 9
	return 0.0;

	#elif init_flag == 10
	if (x<1.0/6.0 + 0.57735026919*y) return -8.25*0.5;
	else                             return 0.0;

	#elif init_flag == 11
	return 0.0;


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
	return 0.0;

	#elif init_flag == 6
	return 0.0;

	#elif init_flag == 7
	return 0.0;

	#elif init_flag == 8
	return 0.0;

	#elif init_flag == 9
	return 0.0;

	#elif init_flag == 10
	return 0.0;

	#elif init_flag == 11
	return 0.0;


	#endif
}

__host__ __device__ double get_r_dust(double x, double y, double z)
{
	#if init_flag == 2
	return D_G_ratio*pow(x/planet_radius,-p_alpha)*get_vertical(x,y,z);

	#elif init_flag == 3
	return D_G_ratio*pow(x*sin(z)/planet_radius,-p_alpha)*get_vertical(x,y,z);

	#else
	return 1.0;

	#endif
}

__host__ __device__ double get_u_dust(double x, double y, double z)
{
	#if init_flag == 2
	double f_gas = get_v(x, y, z)/sqrt(1.0/x);
	double tmp = Stokes*Stokes*f_gas;

	double L = ((1.0+tmp)/tmp) * (1.0 - sqrt(1.0 - 2.0*tmp*(1.0-f_gas)/(1.0 + tmp)/(1.0 + tmp)));
	double u0 = Stokes*(L*L - 2.0*L);
	double u = (1.0 - sqrt(1.0 - 2.0*Stokes*u0))/Stokes;
	
	L = (1.0 - f_gas + Stokes*u/2.0)/(1.0 + Stokes*u/2.0);
	u0 = Stokes*(L*L - 2.0*L);
	u = (1.0 - sqrt(1.0 - 2.0*Stokes*u0))/Stokes;


	return u*sqrt(1.0/x);

	#elif init_flag == 3
	double f_gas = get_v(x, y, z)/sqrt(1.0/x);
	double tmp = Stokes*Stokes*f_gas;

	double L = ((1.0+tmp)/tmp) * (1.0 - sqrt(1.0 - 2.0*tmp*(1.0-f_gas)/(1.0 + tmp)/(1.0 + tmp)));
	double u0 = Stokes*(L*L - 2.0*L);
	double u = (1.0 - sqrt(1.0 - 2.0*Stokes*u0))/Stokes;
	
	L = (1.0 - f_gas + Stokes*u/2.0)/(1.0 + Stokes*u/2.0);
	u0 = Stokes*(L*L - 2.0*L);
	u = (1.0 - sqrt(1.0 - 2.0*Stokes*u0))/Stokes;
	#else
	return 0.0;

	#endif

}

__host__ __device__ double get_v_dust(double x, double y, double z)
{
	#if init_flag == 2
	return sqrt(1.0/x);

	#elif init_flag == 3
	return sqrt(1.0/x);

	#else
	return 0.0;

	#endif

}


__host__ __device__ double get_w_dust(double x, double y, double z)
{
	#if init_flag == 2
	return 0.0;

	#elif init_flag == 3
	return 0.0;

	#else
	return 0.0;

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

__host__ __device__ Dust init_CD(double x, double y, double z)
{
	Dust Q;

	Q.r = get_r_dust(x,y,z);
	Q.u = get_u_dust(x,y,z);
	Q.v = get_v_dust(x,y,z);
	Q.w = get_w_dust(x,y,z);

	return Q;
}

__host__ __device__ Cell init_C(double x0, double x1, double y0, double y1, double z0, double z1)
{
	Cell L, C, R;
	L = init_C(x0,y0,z0);
	C = init_C(0.5*(x0+x1),0.5*(y0+y1),0.5*(z0+z1));
	R = init_C(x1,y1,z1);
	L.multiply(1.0/6.0);
	C.multiply(4.0/6.0);
	R.multiply(1.0/6.0);

	L.add(C);
	L.add(R);

	return L;
}

void make_grid(double* a, double* v, double amin, double amax, int res, int pad, int geom, int grid)
{
	if      (grid==0) linspace(&a[pad], amin, amax, res+1);
	else if (grid==1) logspace(&a[pad], amin, amax, res+1);
	else if (grid==2) nonuspace(&a[pad], amin, amax, res+1);
	else if (grid==3) nonuspace_half(&a[pad], amin, amax, res+1);
	else if (grid==4) nonuspace_mix(&a[pad], amin, amax, res+1);

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

void fill_grid(Grid G)
{
	srand(0);

	int ind;
	for (int i=0; i<G.xarr; i++)
	for (int j=0; j<G.yarr; j++)
	for (int k=0; k<G.zarr; k++)
	{
		ind = G.get_ind(i,j,k);
		//G.C[ind] = init_C(G.get_xa(i), G.get_xa(i+1), G.get_ya(j), G.get_ya(j+1), G.get_za(k), G.get_za(k+1));
		G.C[ind] = init_C(G.get_xc(i), G.get_yc(j), G.get_zc(k));
		#ifdef dust_flag
		G.CD[ind] = init_CD(G.get_xc(i), G.get_yc(j), G.get_zc(k));
		#endif
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
		fill_grid(dev[n]);
	}

	return;
}
