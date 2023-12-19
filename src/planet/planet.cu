#include "parameters.h"
#include "structs.h"
#include <iostream>
#include <fstream>
#include <sstream>

__host__ __device__ double ramp_function(double t, double tau, double val)
{
	if (t<tau)
	{
		return val*sin(hpi*t/tau);
	}
	else
	{
		return val;
	}
}

__global__ void planet_evo(body* planets, double time, double dt)
{
	int n = threadIdx.x;
	planets[n].m = ramp_function(time, ramp_time, planet_mass);

	double rad;
	double fx, fy, fz;

	if (planet_ecc==0.0)
	{
		rad = planets[n].x;
		planets[n].y = fmod((planets[n].vy/rad/rad-frame_omega)*time+pi,twopi);
	}
	else
	{
		rad = planets[n].x;
		fx = planets[n].vy*planets[n].vy/rad/rad/rad - 1.0/rad/rad;
		fy = 0.0;
		fz = 0.0;

		planets[n].vx += 0.5*fx*dt;
		planets[n].vy += 0.5*fy*dt;
		planets[n].vz += 0.5*fz*dt;
	
		planets[n].x += planets[n].vx*dt;
		planets[n].y += (planets[n].vy/planets[n].x/rad-frame_omega)*dt;
		planets[n].z += planets[n].vz*dt;
	
		rad = planets[n].x;
		fx = planets[n].vy*planets[n].vy/rad/rad/rad - 1.0/rad/rad;
		fy = 0.0;
		fz = 0.0;

		planets[n].vx += 0.5*fx*dt;
		planets[n].vy += 0.5*fy*dt;
		planets[n].vz += 0.5*fz*dt;
	}

	return;
}

__host__ __device__ double ecc_anomaly(double n, double e, double t)
{
	double l = fmod(n*t,twopi);
	double D = -e*sin(l);

	if      (e==0.0)   return l;
	else if (l==0.0)   return l;
	else if (l==pi)    return l;
	else if (l==twopi) return l;

	double tor = 1.0e-12;
	int count = 0;

	double x, dx, d2x, eps, step;
	double esinu, ecosu;

	do
	{
		sincos(l-D, &esinu, &ecosu);
		esinu *= e;
		ecosu *= e;
    
        	x   = D   + esinu;
        	dx  = 1.0 - ecosu;
        	d2x = -esinu;

		eps = 2.0*x*d2x/dx/dx;
        
	        if (fabs(eps)<1.0)
		{
			step = -(x/dx)*(1.0 + eps/4.0 + eps*eps/8.0 + 5.0*eps*eps*eps/64.0);
		}
		else
		{
			step = copysign(fmin(fabs(x/dx),hpi/2.0),-x);
		}
        
		D += step;
		count++;
        
	}while (fabs(step)>tor || count < 10);

	return l-D;
}

__host__ __device__ void orbit_solution(body* planet, double t)
{
	double a = (*planet).a;
	double e = (*planet).e;
	double n = sqrt(1.0/a/a/a);
	double efac = sqrt(1.0-e*e);
	double B = e/(1.0+efac);

	double u, cosu, sinu;
	u = ecc_anomaly(n,e,t);
	sincos(u, &sinu, &cosu);

	double f, ecosf, esinf;
	f = u + 2.0*atan(B*sinu/(1.0-B*cosu));
	sincos(f, &esinf, &ecosf);

	esinf *= e;
	ecosf *= e;

	(*planet).m  = ramp_function(t, ramp_time, planet_mass);
	(*planet).x  = a*efac*efac/(1.0+ecosf);
	(*planet).y  = f-fmod(frame_omega*t,twopi)+pi;
	(*planet).z  = hpi;
    	(*planet).vx = n*a*esinf/efac;
    	(*planet).vy = n*(1.0+ecosf)*(1.0+ecosf)/efac/efac/efac;
	(*planet).vz = 0.0;

	return;
}

__host__ __device__ void two_body_solution(body* planet, double t)
{
	double a = planet[0].a;
	double e = planet[0].e;
	double n = sqrt(1.0/a/a/a);
	double efac = sqrt(1.0-e*e);
	double B = e/(1.0+efac);

	double u, cosu, sinu;
	u = ecc_anomaly(n,e,t);
	sincos(u, &sinu, &cosu);

	double f, ecosf, esinf;
	f = u + 2.0*atan(B*sinu/(1.0-B*cosu));
	sincos(f, &esinf, &ecosf);

	esinf *= e;
	ecosf *= e;

	double reduced_r = a*efac*efac/(1.0+ecosf);
	double reduced_p = f-fmod(frame_omega*t,twopi)+pi;
	double mp = ramp_function(t, ramp_time, planet_mass);
	
	planet[0].m  = 1.0/ (1.0+mp);
	planet[0].x  = reduced_r * mp / (1.0+mp);
	planet[0].y  = fmod(reduced_p+pi,twopi);
	planet[0].z  = hpi;
    	planet[0].vx = n*a*esinf/efac * mp / (1.0+mp);
    	planet[0].vy = n*(1.0+ecosf)*(1.0+ecosf)/efac/efac/efac;
	planet[0].vz = 0.0;

	planet[1].m  = mp/ (1.0+mp);
	planet[1].x  = reduced_r * 1.0 / (1.0+mp);
	planet[1].y  = reduced_p;
	planet[1].z  = hpi;
    	planet[1].vx = n*a*esinf/efac * 1.0 / (1.0+mp);
    	planet[1].vy = n*(1.0+ecosf)*(1.0+ecosf)/efac/efac/efac;
	planet[1].vz = 0.0;

	return;
}

__global__ void planet_ana(body* planets, double t)
{
	#if twobd_flag == 1
	two_body_solution(planets, t);
	#else
	int i = threadIdx.x;
	orbit_solution(&planets[i],t);
	#endif

	return;
}

void evolve_planet(Grid* dev, double time, double dt)
{
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);
		//planet_evo<<< 1, n_planet, 0, dev[n].stream >>> (dev[n].planets, time, dt);
		#if twobd_flag == 1
		planet_ana<<< 1, 1, 0, dev[n].stream >>> (dev[n].planets, time);
		#else
		planet_ana<<< 1, n_planet, 0, dev[n].stream >>> (dev[n].planets, time);
		#endif
	}
	return;
}

double get_rs(double a, double mp)
{
	double rs;
	#if ndim==2
	rs = fmin(sc_h,fmax(sc_h/2.0,pow(planet_mass/3.0,1.0/3.0)/4.0))*a;
	#else
	rs = (pow(planet_mass/3.0,1.0/3.0)/4.0)*a;
	#endif
	//rs = 0.0003*a;
	return rs;
}

void init_planet(Grid* G, double time)
{
	printf("planet initialization begins...\n");

	#if twobd_flag == 1
	double e = planet_ecc;
	double a = planet_radius;
	for (int i=0; i<ndev; i++)
	{	
		G[i].planets[0].a = a;
		G[i].planets[0].e = e;

		G[i].planets[1].a = a;
		G[i].planets[1].e = e;

		two_body_solution(G[i].planets, time);

		G[i].planets[0].rs = 0.0;
		G[i].planets[1].rs = get_rs(a, planet_mass);
	}
	#else
	double e = planet_ecc;
	double a;
	for (int i=0; i<ndev; i++)
	{
		for (int n=0; n<n_planet; n++)
		{
			e = planet_ecc;
			if (n==0) a = planet_radius;
			else if (n==n_planet-1) a = 5.0*planet_radius;
			else
			{
				if (n_planet==3) a = 2.236068*planet_radius;
				else if (n_planet==4)
				{
					if (n==1) a = 1.71*planet_radius;
					if (n==2) a = 2.924*planet_radius;
				}
				else if (n_planet==5)
				{
					if (n==1) a = 1.495349*planet_radius;
					if (n==2) a = 2.2360685*planet_radius;
					if (n==3) a = 3.3437*planet_radius;
				}
				else if (n_planet==6)
				{
					if (n==1) a = 1.37973*planet_radius;
					if (n==2) a = 1.903654*planet_radius;
					if (n==3) a = 2.62652784*planet_radius;
					if (n==4) a = 3.6239*planet_radius;
				}
				else if (n_planet==7)
				{
					if (n==1) a = 1.30766044*planet_radius;
					if (n==2) a = 1.70997585*planet_radius;
					if (n==3) a = 2.23606792*planet_radius;
					if (n==4) a = 2.92401782*planet_radius;
					if (n==5) a = 3.82362215*planet_radius;
				}
			}
			
			G[i].planets[n].a = a;
			G[i].planets[n].e = e;

			orbit_solution(&(G[i].planets[n]), time);

			G[i].planets[n].rs = get_rs(a,planet_mass);
/*
			G[i].planets[n].m = ramp_function(time, ramp_time, planet_mass);
			G[i].planets[n].x = a*(1.0-e);
			G[i].planets[n].y = pi;
			G[i].planets[n].z = hpi;

			G[i].planets[n].vx = 0.0;
			G[i].planets[n].vy = sqrt((1.0+e)/(1.0-e)/a)*G[i].planets[n].x;
			G[i].planets[n].vz = 0.0;
*/
		}
	}
	#endif

	printf("planet initialized.\n");

	return;
}
