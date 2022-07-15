__device__ double get_tc(double rad)
{
	double tc;
	tc = pow(rad,1.5)*beta_cool;

	return tc;
}

__device__ double cool(double x, double y, double z, double p, double rho, double dt)
{
	#if geomx==0
	double rad = 1.0;
	#elif geomx==1
	double rad = x;
	#elif geomx==2
	double rad = x*sin(z);
	#endif

	double tc = get_tc(rad);
	double p0 = get_cs2(x,y,z)*rho;
	return (p0-p)*(1.0-exp(-dt/tc));	
}

__global__ void Newtonian_Cooling(Grid G, Cell* in, Cell* out, double dt)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;

	double tc, r, rad_cyl;
	double e, e0;
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

		#if geomx == 0
		rad_cyl = 1.0;
		#elif geomx == 1
		rad_cyl = xc;
		#elif geomx == 2
		rad_cyl = xc * sin(zc);
		#endif

		r = in[ind].r;
		tc = get_tc(rad_cyl);

		e  = in[ind].p/r;
		e0 = get_cs2(xc,yc,zc);

		out[ind].p = r*e0 + r*(e-e0)*exp(-dt/tc);
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
