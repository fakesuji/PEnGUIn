__device__ double get_tc(double rad)
{
	double tc;
	tc = pow(rad,1.5)*beta_cool;

	return tc;
}

__device__ double get_p_cool(double xc, double yc, double zc, double p, double r, double dt)
{
	double rad_cyl;

	#if geomx == 0
	rad_cyl = 1.0;
	#elif geomx == 1
	rad_cyl = xc;
	#elif geomx == 2
	rad_cyl = xc * sin(zc);
	#endif

	double tc = get_tc(rad_cyl);

	double e  = p/r;
	double e0 = get_cs2(xc,yc,zc);

	return r*(e0 + (e-e0)*exp(-dt/tc));
}

__global__ void Newtonian_Cooling(Grid G, Cell* in, Cell* out, double dt)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;

	double xc, yc, zc;

	int ind;

	if (i>=xpad && i<G.xarr-xpad)
	if (j>=ypad && j<G.yarr-ypad)
	if (k>=zpad && k<G.zarr-zpad)
	{		
		ind = G.get_ind(i,j,k);

		xc = G.get_xc(i);
		yc = G.get_yc(j);
		zc = G.get_zc(k);

		out[ind].p = get_p_cool(xc, yc, zc, in[ind].p, in[ind].r, dt);
	}

	return;
}

void cooling(Grid* dev, double dt)
{
	int nx, ny, nz;
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		Newtonian_Cooling<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 0, dev[n].stream >>> (dev[n],dev[n].C,dev[n].C,dt);
	}
	return;
}
