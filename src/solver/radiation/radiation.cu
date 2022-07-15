__global__ void get_tau(Grid G, double kappa)
{
	int i;
	int j = threadIdx.x + blockIdx.x*blockDim.x + ypad;
	int k = threadIdx.y + blockIdx.y*blockDim.y + zpad;

	int ix;

	double tmp = 0.0;

	for (i=0; i<G.xarr-1; i++)
	{
		ix = G.get_ind(i,j,k);

		G.ext[ix]  = tmp;
		tmp = kappa * G.C[ix].r * 0.5*(G.get_xa(i+1)-G.get_xa(i));
		G.ext[ix] += tmp;

		G.ext[G.get_ind(i+1,j,k)] += G.ext[ix]; 
	}

	i = G.xarr-1;
	ix = G.get_ind(i,j,k);

	G.ext[ix]  = tmp;
	tmp = kappa * G.C[ix].r * 0.5*(G.get_xa(i+1)-G.get_xa(i));
	G.ext[ix] += tmp;

	return;
}

void compute_extinction(Grid* dev, double kappa)
{
	int ny, nz;
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);
		ny = dev[n].yres;
		nz = dev[n].zres;

		get_tau<<< dim3(ny/x_ydiv,nz/x_zdiv), dim3(x_ydiv,x_zdiv), 0, dev[n].stream >>> (dev[n],kappa);
	}
	return;
}
