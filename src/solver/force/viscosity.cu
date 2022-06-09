__global__ void viscous_tensor(Grid G, Cell* C)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + 1;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;
	#if ndim>1
	j++;
	#endif
	#if ndim>2
	k++;
	#endif
	if (i>=G.xarr-1) return;

	int s;
	#if ndim==1
	s = 1;
	#elif ndim==2
	s = 3;
	#elif ndim==3
	s = 6;
	#endif

	int ind;
	double up, vp, wp;
	double um, vm, wm;

	double r = G.get_xc(i);
	double dr = G.get_xc(i+1)-G.get_xc(i-1);

	ind = G.get_ind(i+1, j, k);
	up = C[ind].u;
	vp = C[ind].v;
	wp = C[ind].w;

	ind = G.get_ind(i-1, j, k);
	um = C[ind].u;
	vm = C[ind].v;
	wm = C[ind].w;

	double dudx = (up-um)/dr;
	#if ndim>1
	double dvdx = (vp-vm)/dr;
	#endif
	#if ndim>2
	double dwdx = (wp-wm)/dr;
	#endif

	#if ndim>1
	double dphi = G.get_yc(j+1)-G.get_yc(j-1);

	ind = G.get_ind(i, j+1, k);
	up = C[ind].u;
	vp = C[ind].v;
	wp = C[ind].w;

	ind = G.get_ind(i, j-1, k);
	um = C[ind].u;
	vm = C[ind].v;
	wm = C[ind].w;

	double dudy = (up-um)/dphi/r;
	double dvdy = (vp-vm)/dphi/r;
	#if ndim>2
	double dwdy = (wp-wm)/dphi/r;
	#endif
	#endif

	#if ndim>2
	double dz   = G.get_zc(k+1)-G.get_zc(k-1);

	ind = G.get_ind(i, j, k+1);
	up = C[ind].u;
	vp = C[ind].v;
	wp = C[ind].w;

	ind = G.get_ind(i, j, k-1);
	um = C[ind].u;
	vm = C[ind].v;
	wm = C[ind].w;

	double dudz = (up-um)/dz/r;
	double dvdz = (vp-vm)/dz/r;
	double dwdz = (wp-wm)/dz/r;
	#endif

	ind = G.get_ind(i, j, k);
	double rho = C[ind].r;
	double u = C[ind].u;
	double v = C[ind].v;

	double nr = get_nu(r,G.get_yc(j),G.get_zc(k))*rho;

	G.vis_tensor[s*ind+0] = nr*(dudx);

	#if ndim>1
	G.vis_tensor[s*ind+1] = nr*(dudy + dvdx - v/r);
	G.vis_tensor[s*ind+2] = nr*(dvdy + u/r);
	#endif

	#if ndim>2
	G.vis_tensor[s*ind+3] = nr*(dudz + dwdx);
	G.vis_tensor[s*ind+4] = nr*(dvdz + dwdy);
	G.vis_tensor[s*ind+5] = nr*(dwdz);
	#endif

	//printf("%f,%e,%e\n",r,nr,dudx);
	return;
}

__device__ double viscous_fx(Grid G, Cell* C, int i, int j, int k)
{
	int s;
	#if ndim==1
	s = 1;
	#elif ndim==2
	s = 3;
	#elif ndim==3
	s = 6;
	#endif

	int ind;
	double fx;

	ind = G.get_ind(i+1, j, k);
	double* tf = &G.vis_tensor[s*ind];

	ind = G.get_ind(i-1, j, k);
	double* tb = &G.vis_tensor[s*ind];

	#if ndim>1
	ind = G.get_ind(i, j+1, k);
	double* tr = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j-1, k);
	double* tl = &G.vis_tensor[s*ind];
	#endif

	#if ndim>2
	ind = G.get_ind(i, j, k+1);
	double* tt = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j, k-1);
	double* td = &G.vis_tensor[s*ind];
	#endif

	ind = G.get_ind(i, j, k);
	double* t0 = &G.vis_tensor[s*ind];

	double r = G.get_xc(i);
	double dr   = G.get_xc(i+1)-G.get_xc(i-1);
	fx = 2.0*t0[0]/r + 2.0*(tf[0]-tb[0])/dr;

	#if ndim>1
	double dphi = G.get_yc(j+1)-G.get_yc(j-1);
	fx += (tr[1]-tl[1])/dphi - 2.0*t0[2]/r;
	#endif

	#if ndim>2
	double dz   = G.get_zc(k+1)-G.get_zc(k-1);
	fx += (tt[3]-td[3])/dz;
	#endif

	return fx;
}

__device__ double viscous_fy(Grid G, Cell* C, int i, int j, int k)
{
	int s;
	#if ndim==1
	s = 1;
	#elif ndim==2
	s = 3;
	#elif ndim==3
	s = 6;
	#endif

	int ind;
	double fy;

	ind = G.get_ind(i+1, j, k);
	double* tf = &G.vis_tensor[s*ind];

	ind = G.get_ind(i-1, j, k);
	double* tb = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j+1, k);
	double* tr = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j-1, k);
	double* tl = &G.vis_tensor[s*ind];

	#if ndim>2
	ind = G.get_ind(i, j, k+1);
	double* tt = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j, k-1);
	double* td = &G.vis_tensor[s*ind];
	#endif

	ind = G.get_ind(i, j, k);
	double* t0 = &G.vis_tensor[s*ind];

	double r = G.get_xc(i);
	double dr   = G.get_xc(i+1)-G.get_xc(i-1);
	double dphi = G.get_yc(j+1)-G.get_yc(j-1);

	fy = 2.0*t0[1]/r + (tf[1]-tb[1])/dr + 2.0*(tr[2]-tl[2])/dphi;

	#if ndim>2
	double dz   = G.get_zc(k+1)-G.get_zc(k-1);
	fy += (tt[4]-td[4])/dz;
	#endif

	return fy;
}

__device__ double viscous_fz(Grid G, Cell* C, int i, int j, int k)
{
	int s;
	#if ndim==1
	s = 1;
	#elif ndim==2
	s = 3;
	#elif ndim==3
	s = 6;
	#endif

	int ind;
	double fz;

	ind = G.get_ind(i+1, j, k);
	double* tf = &G.vis_tensor[s*ind];

	ind = G.get_ind(i-1, j, k);
	double* tb = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j+1, k);
	double* tr = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j-1, k);
	double* tl = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j, k+1);
	double* tt = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j, k-1);
	double* td = &G.vis_tensor[s*ind];

	ind = G.get_ind(i, j, k);
	double* t0 = &G.vis_tensor[s*ind];

	double r = G.get_xc(i);
	double dr   = G.get_xc(i+1)-G.get_xc(i-1);
	double dphi = G.get_yc(j+1)-G.get_yc(j-1);
	double dz   = G.get_zc(k+1)-G.get_zc(k-1);

	fz = t0[3]/r + (tf[3]-tb[3])/dr + (tr[4]-tl[4])/dphi + 2.0*(tt[5]-td[5])/dz;

	return fz;
}

__device__ double viscous_heat(Grid G, Cell* C, int i, int j, int k)
{
	int s;
	#if ndim==1
	s = 1;
	#elif ndim==2
	s = 3;
	#elif ndim==3
	s = 6;
	#endif

	int ind;
	double H;

	ind = G.get_ind(i, j, k);
	double* T = &G.vis_tensor[s*ind];

	double rho = C[ind].r;
	double nr = get_nu(G.get_xc(i),G.get_yc(j),G.get_zc(k))*rho;

	H = 2.0*T[0]*T[0] + 2.0*T[2]*T[2] + 2.0*T[5]*T[5] + T[1]*T[1] + T[3]*T[3] + T[4]*T[4];
	H /= nr;

	return H;
}

void viscosity_tensor_evaluation1(Grid* dev)
{
	int nx, ny, nz;

	#if ndim>2
	boundz(dev);
	#endif
	#if ndim>1
	boundy(dev);
	#endif
	boundx(dev);
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xarr-2;
		#if ndim>1
		ny = dev[n].yarr-2;
		#else
		ny = 1;
		#endif
		#if ndim>2
		nz = dev[n].zarr-2;
		#else
		nz = 1;
		#endif
		viscous_tensor<<< dim3((nx+255)/256,ny,nz), 256, 0, dev[n].stream >>>(dev[n], dev[n].C);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
}

void viscosity_tensor_evaluation2(Grid* dev)
{
	int nx, ny, nz;

	#if ndim>2
	boundz2(dev);
	#endif
	#if ndim>1
	boundy2(dev);
	#endif
	boundx2(dev);
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xarr-2;
		#if ndim>1
		ny = dev[n].yarr-2;
		#else
		ny = 1;
		#endif
		#if ndim>2
		nz = dev[n].zarr-2;
		#else
		nz = 1;
		#endif
		viscous_tensor<<< dim3((nx+255)/256,ny,nz), 256, 0, dev[n].stream >>>(dev[n], dev[n].T);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
}

