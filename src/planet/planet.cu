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
	#if geomx==0
	rad = 1.0;
	#else
	rad = planets[n].x + 0.5*planets[n].vx*dt;
	#endif

	planets[n].x += planets[n].vx*dt;
	planets[n].y += (planets[n].vy-frame_omega*rad)*dt;
	planets[n].z += planets[n].vz*dt;
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
	for (int i=0; i<ndev; i++)
	{
		for (int n=0; n<n_planet; n++)
		{
			G[i].planets[n].m = ramp_function(time, ramp_time, planet_mass);
			G[i].planets[n].x = planet_radius;
			G[i].planets[n].y = pi;
			G[i].planets[n].z = hpi;

			G[i].planets[n].vx = 0.0;
			G[i].planets[n].vy = pow(G[i].planets[n].x,-1.5);
			G[i].planets[n].vz = 0.0;
			#if ndim==2
			G[i].planets[n].rs = 0.5*sc_h;
			#else
			G[i].planets[n].rs = 0.5*pow(planet_mass,1.0/3.0);
			#endif
		}
	}
	return;
}
