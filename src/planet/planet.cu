#include "parameters.h"
#include "structs.h"

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
	//	planets[n].y += (planets[n].vy/rad/rad-frame_omega)*dt;
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

void evolve_planet(Grid* dev, double time, double dt)
{
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);
		planet_evo<<< 1, n_planet, 0, dev[n].stream >>> (dev[n].planets, time+dt, dt);
	}
	return;
}

void init_planet(Grid* G, double time)
{
	double e = planet_ecc;
	double a;
	for (int i=0; i<ndev; i++)
	{
		for (int n=0; n<n_planet; n++)
		{
			e = planet_ecc;
			if (n==0) a = planet_radius;//planet_radius;
			else if (n==n_planet-1) a = 5.0;
			else
			{
				if (n_planet==3) a = 2.236068;
				else if (n_planet==4)
				{
					if (n==1) a = 1.71;
					if (n==2) a = 2.924;
				}
				else if (n_planet==5)
				{
					if (n==1) a = 1.495349;
					if (n==2) a = 2.2360685;
					if (n==3) a = 3.3437;
				}
				else if (n_planet==6)
				{
					if (n==1) a = 1.37973;
					if (n==2) a = 1.903654;
					if (n==3) a = 2.62652784;
					if (n==4) a = 3.6239;
				}
				else if (n_planet==7)
				{
					if (n==1) a = 1.30766044;
					if (n==2) a = 1.70997585;
					if (n==3) a = 2.23606792;
					if (n==4) a = 2.92401782;
					if (n==5) a = 3.82362215;
				}
			}
			
			G[i].planets[n].m = ramp_function(time, ramp_time, planet_mass);
			G[i].planets[n].x = a*(1.0-e);
			G[i].planets[n].y = pi;
			G[i].planets[n].z = hpi;

			G[i].planets[n].vx = 0.0;
			G[i].planets[n].vy = sqrt((1.0+e)/(1.0-e)/a)*G[i].planets[n].x;
			G[i].planets[n].vz = 0.0;
			#if ndim==2
			G[i].planets[n].rs = fmin(sc_h,fmax(sc_h/2.0,pow(planet_mass/3.0,1.0/3.0)/4.0))*a;
			#else
			G[i].planets[n].rs = (pow(planet_mass/3.0,1.0/3.0)/4.0)*a;
			#endif
		}
	}
	return;
}
