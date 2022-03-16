__device__ double ramp_function(double t, double tau, double val)
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
