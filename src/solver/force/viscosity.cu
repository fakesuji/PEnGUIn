__global__ void viscous_tensor(Grid G)
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
	#elif nidm==3
	s = 6;
	#endif

	int ind;
	double up, vp, wp;
	double um, vm, wm;

	double r = G.get_xc(i);
	double dr = G.get_xc(i+1)-G.get_xc(i-1);

	ind = G.get_ind(i+1, j, k);
	up = G.C[ind].u;
	vp = G.C[ind].v;
	wp = G.C[ind].w;

	ind = G.get_ind(i-1, j, k);
	um = G.C[ind].u;
	vm = G.C[ind].v;
	wm = G.C[ind].w;

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
	up = G.C[ind].u;
	vp = G.C[ind].v;
	wp = G.C[ind].w;

	ind = G.get_ind(i, j-1, k);
	um = G.C[ind].u;
	vm = G.C[ind].v;
	wm = G.C[ind].w;

	double dudy = (up-um)/dphi/r;
	double dvdy = (vp-vm)/dphi/r;
	#if ndim>2
	double dwdy = (wp-wm)/dphi/r;
	#endif
	#endif

	#if ndim>2
	double dz   = G.get_zc(k+1)-G.get_zc(k-1);

	ind = G.get_ind(i, j, k+1);
	up = G.C[ind].u;
	vp = G.C[ind].v;
	wp = G.C[ind].w;

	ind = G.get_ind(i, j, k-1);
	um = G.C[ind].u;
	vm = G.C[ind].v;
	wm = G.C[ind].w;

	double dudz = (up-um)/dz/r;
	double dvdz = (vp-vm)/dz/r;
	double dwdz = (wp-wm)/dz/r;
	#endif

	ind = G.get_ind(i, j, k);
	double rho = G.C[ind].r;
	double u = G.C[ind].u;
	double v = G.C[ind].v;
	//double w = G.C[ind].w;

	double nr = get_nu(r,G.get_yc(j),G.get_zc(k))*rho;

	G.vis_tensor[s*ind+0] = 2.0*nr*(dudx);
	#if ndim>1
	G.vis_tensor[s*ind+1] =     nr*(dudy + dvdx - v/r);
	G.vis_tensor[s*ind+2] = 2.0*nr*(dvdy + u/r);
	#endif
	#if ndim>2
	G.vis_tensor[s*ind+3] =     nr*(dudz + dwdx);
	G.vis_tensor[s*ind+4] =     nr*(dvdz + dwdy);
	G.vis_tensor[s*ind+5] = 2.0*nr*(dwdz);
	#endif

	//printf("%f,%e,%e\n",r,nr,dudx);
	return;
}

__global__ void viscous_force(Grid G, double dt)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;

	int s;
	#if ndim==1
	s = 1;
	#elif ndim==2
	s = 3;
	#elif nidm==3
	s = 6;
	#endif

	int ind;
	double fx, fy, fz;

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
	double rho = G.C[ind].r;
	double dr   = G.get_xc(i+1)-G.get_xc(i-1);
	fx = t0[0]/r + (tf[0]-tb[0])/dr;
	fy = 0.0;
	fz = 0.0;

	#if ndim>1
	double dphi = G.get_yc(j+1)-G.get_yc(j-1);
	fx += (tr[1]-tl[1])/dphi - t0[2]/r;
	fy += 2.0*t0[1]/r + (tf[1]-tb[1])/dr + (tr[2]-tl[2])/dphi;
	#endif

	#if ndim>2
	double dz   = G.get_zc(k+1)-G.get_zc(k-1);
	fx += (tt[3]-td[3])/dz;
	fy += (tt[4]-td[4])/dz;
	fz += t0[3]/r + (tf[3]-tb[3])/dr + (tr[4]-tl[4])/dphi + (tt[5]-td[5])/dz;
	#endif

	//if(j==ypad) printf("%f,%e,%e\n",r,fx,fy);

	G.C[ind].u += fx*dt/rho;
	G.C[ind].v += fy*dt/rho;
	G.C[ind].w += fz*dt/rho;

	return;
}

void viscosity_tensor_evaluation(Grid* dev)
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
		//printf("tensor calculation\n");
		viscous_tensor<<< dim3((nx+255)/256,ny,nz), 256, 0, dev[n].stream >>>(dev[n]);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
}

void apply_viscosity(Grid* dev, double dt)
{
	int nx, ny, nz;

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres/x_xdiv;
		ny = dev[n].yres/x_ydiv;
		nz = dev[n].zres/x_zdiv;
		viscous_force <<< dim3(nx,ny,nz), dim3(x_xdiv,x_ydiv,x_zdiv), 0, dev[n].stream >>>(dev[n],dt);
	}
	for (int n=0; n<ndev; n++) cudaStreamSynchronize(dev[n].stream);
}
