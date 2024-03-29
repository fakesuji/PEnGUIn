__device__ Cell advection(int geom, double* xa, double* dx, double* dv, int npad, double* r, double* e, double* u, double* v, double* w, double speed, double dt)
{
	Cell Del;

	int imax = blockDim.x;
	int i = threadIdx.x;

	double xl, xr, x0;
	double q, ql, qr;
	double par[4];

	#if recon_flag==6
	extern __shared__ double share[];
	double* flat = &share[0];
	int is = i + imax*threadIdx.y;
	flat[is] = 0.0;
	__syncthreads();
	#endif

	if (i>=npad && i<imax+1-npad)
	{
		if (speed>0.0)
		{
			xl = xa[i-1];
			xr = xa[i];
			x0 = xr - speed*dt;
			dimensionless_x(xl,x0,xr,q,ql,qr);
	
			get_CON_parameters(i-1, geom, xa, dx, dv, r, par);
			Del.r = get_CON_aveR(geom, q, par, ql, qr);
	
			get_CON_parameters(i-1, geom, xa, dx, dv, e, par);
			Del.p = get_CON_aveR(geom, q, par, ql, qr);

			get_PRM_parameters(i-1, geom, xa, dx, dv, u, par);
			Del.u = get_PRM_aveR(geom, q, par, ql, qr)*Del.r;

			get_PRM_parameters(i-1, geom, xa, dx, dv, v, par);
			Del.v = get_PRM_aveR(geom, q, par, ql, qr)*Del.r;

			get_PRM_parameters(i-1, geom, xa, dx, dv, w, par);
			Del.w = get_PRM_aveR(geom, q, par, ql, qr)*Del.r;
		}
		else
		{
			xl = xa[i];
			xr = xa[i+1];
			x0 = xl - speed*dt;
			dimensionless_x(xl,x0,xr,q,ql,qr);
	
			get_CON_parameters(i, geom, xa, dx, dv, r, par);
			Del.r = get_CON_aveL(geom, q, par, ql, qr);
	
			get_CON_parameters(i, geom, xa, dx, dv, e, par);
			Del.p = get_CON_aveL(geom, q, par, ql, qr);
	
			get_PRM_parameters(i, geom, xa, dx, dv, u, par);
			Del.u = get_PRM_aveL(geom, q, par, ql, qr)*Del.r;

			get_PRM_parameters(i, geom, xa, dx, dv, v, par);
			Del.v = get_PRM_aveL(geom, q, par, ql, qr)*Del.r;

			get_PRM_parameters(i, geom, xa, dx, dv, w, par);
			Del.w = get_PRM_aveL(geom, q, par, ql, qr)*Del.r;
		}
	}
	Del.multiply(speed);
	__syncthreads();

	extern __shared__ double share[];
	double* tmp1 = &share[0];
	double gfac = eval_gfac(geom,xa[i]);

	Del.r = net_flux(tmp1,Del.r*gfac);
	__syncthreads();
	Del.p = net_flux(tmp1,Del.p*gfac);
	__syncthreads();
	Del.u = net_flux(tmp1,Del.u*gfac);
	__syncthreads();
	Del.v = net_flux(tmp1,Del.v*gfac);
	__syncthreads();
	Del.w = net_flux(tmp1,Del.w*gfac);
	__syncthreads();

	return Del;
}

__device__ Cell update_Cell_advect(double x, double y, double z, Cell Q, Cell D)
{
	Q.u *= Q.r;
	Q.v *= Q.r;
	Q.w *= Q.r;

	Q.add(D);
	
	Q.u /= Q.r;
	Q.v /= Q.r;
	Q.w /= Q.r;

	#if EOS_flag == 0
	Q.p = get_cs2(x, y, z)*Q.r;
	#endif

	return Q;
}

__global__ void advecty(Grid G, Cell* in, Cell* out, double dt)
{
	__shared__ double ya[y_ythd+1], dy[y_ythd], yv[y_ythd];
	__shared__ double r[y_ythd*y_xdiv], p[y_ythd*y_xdiv], u[y_ythd*y_xdiv], v[y_ythd*y_xdiv], w[y_ythd*y_xdiv];

	int i = threadIdx.x;
	int idy = i + blockIdx.x*y_ydiv;

	int j = threadIdx.y;
	int idx = j + blockIdx.y*y_xdiv + xpad;

	int k = threadIdx.z;
	int idz = k + blockIdx.z*y_zdiv + zpad;

	int ind = idx + G.xarr*(idy + G.yarr*idz);

	if (j==0)
	{
		ya[i] = G.get_ya(idy);
		if (i==blockDim.x-1) ya[i+1] = G.get_ya(idy+1);
		yv[i] = G.get_yv(idy);
	}
	__syncthreads();

	if (j==0) dy[i] = ya[i+1] - ya[i];

	double rad;	
	#if geomy == 3
	rad = G.get_xc(idx);
	#elif geomy == 4
	rad = G.get_xc(idx);
	rad *= sin(G.get_zc(idz));
	#else
	rad = 1.0;
	#endif

	r[i+y_ythd*j] = in[ind].r;
	p[i+y_ythd*j] = in[ind].p;
	u[i+y_ythd*j] = in[ind].v;
	v[i+y_ythd*j] = in[ind].w;
	w[i+y_ythd*j] = in[ind].u;
	__syncthreads();

	/////////////////////////////////////////////////////

	Cell Del;
	Del = advection(geomy, ya, dy, yv, ypad, &r[y_ythd*j], &p[y_ythd*j], &u[y_ythd*j], &v[y_ythd*j], &w[y_ythd*j], G.get_res(idx,idz), dt/rad);
	Del.multiply(dt/yv[i]);

	#if geomy > 2
	Del.multiply(1.0/rad);
	#endif

	double tmp;

	if (i>=ypad && i<y_ythd-ypad)
	{
		tmp = Del.v;
		Del.v = Del.u;
		Del.u = Del.w;
		Del.w = tmp;

		out[ind] = update_Cell_advect(G.get_xc(idx), G.get_yc(idy), G.get_zc(idz), in[ind], Del);
	}

	return;
}

__device__ Dust advection(int geom, double* xa, double* dx, double* dv, int npad, double* r, double* u, double* v, double* w, double speed, double dt)
{
	Dust Del;

	int imax = blockDim.x;
	int i = threadIdx.x;

	double xl, xr, x0;
	double q, ql, qr;
	double par[4];

	#if recon_flag==2
	extern __shared__ double share[];
	double* flat = &share[0];
	int is = i + imax*threadIdx.y;
	flat[is] = 0.0;
	__syncthreads();
	#endif

	if (i>=npad && i<imax+1-npad)
	{
		if (speed>0.0)
		{
			xl = xa[i-1];
			xr = xa[i];
			x0 = xr - speed*dt;
			dimensionless_x(xl,x0,xr,q,ql,qr);
	
			get_CON_parameters(i-1, geom, xa, dx, dv, r, par);
			Del.r = get_CON_aveR(geom, q, par, ql, qr);

			get_PRM_parameters(i-1, geom, xa, dx, dv, u, par);
			Del.u = get_PRM_aveR(geom, q, par, ql, qr)*Del.r;

			get_PRM_parameters(i-1, geom, xa, dx, dv, v, par);
			Del.v = get_PRM_aveR(geom, q, par, ql, qr)*Del.r;

			get_PRM_parameters(i-1, geom, xa, dx, dv, w, par);
			Del.w = get_PRM_aveR(geom, q, par, ql, qr)*Del.r;
		}
		else
		{
			xl = xa[i];
			xr = xa[i+1];
			x0 = xl - speed*dt;
			dimensionless_x(xl,x0,xr,q,ql,qr);
	
			get_CON_parameters(i, geom, xa, dx, dv, r, par);
			Del.r = get_CON_aveL(geom, q, par, ql, qr);
	
			get_PRM_parameters(i, geom, xa, dx, dv, u, par);
			Del.u = get_PRM_aveL(geom, q, par, ql, qr)*Del.r;

			get_PRM_parameters(i, geom, xa, dx, dv, v, par);
			Del.v = get_PRM_aveL(geom, q, par, ql, qr)*Del.r;

			get_PRM_parameters(i, geom, xa, dx, dv, w, par);
			Del.w = get_PRM_aveL(geom, q, par, ql, qr)*Del.r;
		}
	}
	Del.multiply(speed);
	__syncthreads();

	extern __shared__ double share[];
	double* tmp1 = &share[0];
	double gfac = eval_gfac(geom,xa[i]);

	Del.r = net_flux(tmp1,Del.r*gfac);
	__syncthreads();
	Del.u = net_flux(tmp1,Del.u*gfac);
	__syncthreads();
	Del.v = net_flux(tmp1,Del.v*gfac);
	__syncthreads();
	Del.w = net_flux(tmp1,Del.w*gfac);
	__syncthreads();

	return Del;
}

__global__ void advecty(Grid G, Dust* in, Dust* out, double dt)
{
	__shared__ double ya[y_ythd+1], dy[y_ythd], yv[y_ythd];
	__shared__ double r[y_ythd*y_xdiv], u[y_ythd*y_xdiv], v[y_ythd*y_xdiv], w[y_ythd*y_xdiv];

	int i = threadIdx.x;
	int idy = i + blockIdx.x*y_ydiv;

	int j = threadIdx.y;
	int idx = j + blockIdx.y*y_xdiv + xpad;

	int k = threadIdx.z;
	int idz = k + blockIdx.z*y_zdiv + zpad;

	int ind = idx + G.xarr*(idy + G.yarr*idz);

	if (j==0)
	{
		ya[i] = G.get_ya(idy);
		if (i==blockDim.x-1) ya[i+1] = G.get_ya(idy+1);
		yv[i] = G.get_yv(idy);
	}
	__syncthreads();

	if (j==0) dy[i] = ya[i+1] - ya[i];

	double rad;	
	#if geomy == 3
	rad = G.get_xc(idx);
	#elif geomy == 4
	rad = G.get_xc(idx);
	rad *= sin(G.get_zc(idz));
	#else
	rad = 1.0;
	#endif

	r[i+y_ythd*j] = in[ind].r;
	u[i+y_ythd*j] = in[ind].v;
	v[i+y_ythd*j] = in[ind].w;
	w[i+y_ythd*j] = in[ind].u;
	__syncthreads();

	/////////////////////////////////////////////////////

	Dust Del;
	Del = advection(geomy, ya, dy, yv, ypad, &r[y_ythd*j], &u[y_ythd*j], &v[y_ythd*j], &w[y_ythd*j], G.get_res(idx,idz), dt/rad);
	Del.multiply(dt/yv[i]);

	#if geomy > 2
	Del.multiply(1.0/rad);
	#endif

	double tmp;

	if (i>=ypad && i<y_ythd-ypad)
	{
		tmp = Del.v;
		Del.v = Del.u;
		Del.u = Del.w;
		Del.w = tmp;

		out[ind] = update_Dust(G.get_xc(idx), G.get_yc(idy), G.get_zc(idz), in[ind], Del);
	}

	return;
}

void advecty(Grid* dev, double dt)
{
	int nx,ny,nz;

	boundy(dev);
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		advecty<<< dim3(ny/y_ydiv,nx/y_xdiv,nz/y_zdiv), dim3(y_ythd,y_xdiv,y_zdiv), 2*sizeof(double)*y_ythd*y_xdiv*y_zdiv, dev[n].stream >>>
		      (dev[n], dev[n].C, dev[n].T, dt);
		dev[n].CT_change();
	}
}

void advecty_dust(Grid* dev, double dt)
{
	int nx,ny,nz;

	boundy_dust(dev);
	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		advecty<<< dim3(ny/y_ydiv,nx/y_xdiv,nz/y_zdiv), dim3(y_ythd,y_xdiv,y_zdiv), 2*sizeof(double)*y_ythd*y_xdiv*y_zdiv, dev[n].stream >>>
		      (dev[n], dev[n].CD, dev[n].TD, dt);
		dev[n].CT_D_change();
	}
}
