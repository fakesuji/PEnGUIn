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
	for (int i=0; i<ndev; i++)
	{
		for (int n=0; n<n_planet; n++)
		{
			G[i].planets[n].m = ramp_function(time, ramp_time, planet_mass);
			G[i].planets[n].x = planet_radius*(1.0-e);
			G[i].planets[n].y = pi;
			G[i].planets[n].z = hpi;

			G[i].planets[n].vx = 0.0;
			G[i].planets[n].vy = sqrt((1.0+e)/(1.0-e)/planet_radius)*G[i].planets[n].x;
			G[i].planets[n].vz = 0.0;
			#if ndim==2
			G[i].planets[n].rs = fmax(sc_h/2.0,pow(planet_mass,1.0/3.0)/3.0)*planet_radius;
			#else
			G[i].planets[n].rs = (pow(planet_mass,1.0/3.0)/3.0)*planet_radius;
			#endif
		}
	}
	return;
}
