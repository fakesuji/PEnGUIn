__device__ Cell advection(int geom, double* xa, double* dx, double* dv, int npad, double* r, double* e, double* u, double* v, double* w, double speed, double dt)
{
	Cell Del;

	int imax = blockDim.x;
	int i = threadIdx.x;

	double xl, xr, tmp;
	double par[4];

	if (i>=npad && i<imax+1-npad)
	{
		if (speed>0.0)
		{
			xl = xa[i-1];
			xr = xa[i];
			tmp = xr - speed*dt;
	
			get_CON_parameters(i-1, geom, xa, dx, dv, r, par);
			Del.r = get_CON_aveR(geom, xl, tmp, xr, par);
	
			get_CON_parameters(i-1, geom, xa, dx, dv, e, par);
			Del.p = get_CON_aveR(geom, xl, tmp, xr, par);

			get_PRM_parameters(i-1, geom, xa, dx, dv, u, par);
			Del.u = get_PRM_aveR(geom, xl, tmp, xr, par)*Del.r;

			#if ndim>1
			get_PRM_parameters(i-1, geom, xa, dx, dv, v, par);
			Del.v = get_PRM_aveR(geom, xl, tmp, xr, par)*Del.r;
			#else
			Del.v = 0.0;
			#endif

			#if ndim>2
			get_PRM_parameters(i-1, geom, xa, dx, dv, w, par);
			Del.w = get_PRM_aveR(geom, xl, tmp, xr, par)*Del.r;
			#else
			Del.w = 0.0;
			#endif
		}
		else
		{
			xl = xa[i];
			xr = xa[i+1];
			tmp = xl - speed*dt;
	
			get_CON_parameters(i, geom, xa, dx, dv, r, par);
			Del.r = get_CON_aveL(geom, xl, tmp, xr, par);
	
			get_CON_parameters(i, geom, xa, dx, dv, e, par);
			Del.p = get_CON_aveL(geom, xl, tmp, xr, par);
	
			get_PRM_parameters(i, geom, xa, dx, dv, u, par);
			Del.u = get_PRM_aveL(geom, xl, tmp, xr, par)*Del.r;

			#if ndim>1
			get_PRM_parameters(i, geom, xa, dx, dv, v, par);
			Del.v = get_PRM_aveL(geom, xl, tmp, xr, par)*Del.r;
			#else
			Del.v = 0.0;
			#endif

			#if ndim>2
			get_PRM_parameters(i, geom, xa, dx, dv, w, par);
			Del.w = get_PRM_aveL(geom, xl, tmp, xr, par)*Del.r;
			#else
			Del.w = 0.0;
			#endif
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

__global__ void advecy(Grid G, double dt)
{
	__shared__ double ya[y_ythd+1], dy[y_ythd], yv[y_ythd];
	__shared__ double r[y_ythd*y_xdiv], p[y_ythd*y_xdiv], u[y_ythd*y_xdiv], v[y_ythd*y_xdiv], w[y_ythd*y_xdiv];
	#if advection_flag == 1
	__shared__ double s[y_ythd*y_xdiv];
	#endif

	int i = threadIdx.x;
	int idy = i + blockIdx.x*y_ydiv;

	int j = threadIdx.y;
	int idx = j + blockIdx.y*y_xdiv + xpad;

	int k = threadIdx.z;
	int idz = k + blockIdx.z*x_zdiv + zpad;

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

	r[i+y_ythd*j] = G.C[ind].r;
	p[i+y_ythd*j] = G.C[ind].p;
	u[i+y_ythd*j] = G.C[ind].v;
	v[i+y_ythd*j] = G.C[ind].w;
	w[i+y_ythd*j] = G.C[ind].u;
	__syncthreads();

	/////////////////////////////////////////////////////

	Cell Del;
	Del = advection(geomy, ya, dy, yv, ypad, &r[y_ythd*j], &p[y_ythd*j], &u[y_ythd*j], &v[y_ythd*j], &w[y_ythd*j], G.get_res(idx,idz), dt/rad);
	Del.multiply(G.get_xv(idx)*G.get_zv(idz));

	#if geomy > 2
	Del.multiply(1.0/rad);
	#endif

	if (i>=ypad && i<y_ythd-ypad)
	{
		G.F[ind].r = Del.r;
		G.F[ind].p = Del.p;
		G.F[ind].u = Del.w;
		G.F[ind].v = Del.u;
		G.F[ind].w = Del.v;
	}

	return;
}

__global__ void updadv(Grid G, double dt)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x + xpad;
	int j = threadIdx.y + blockIdx.y*blockDim.y + ypad;
	int k = threadIdx.z + blockIdx.z*blockDim.z + zpad;
	double vol;
	Cell Q;
	Cell D;
	int ind;

	if (i>=xpad && i<G.xarr-xpad)
	if (j>=ypad && j<G.yarr-ypad)
	if (k>=zpad && k<G.zarr-zpad)
	{
		ind = i + G.xarr*(j + G.yarr*k);
		vol = G.get_xv(i)*G.get_yv(j)*G.get_zv(k);

		Q.copy(G.C[ind]);
		D.copy(G.F[ind]);		
		D.multiply(dt);

		Q.r *= vol;
		Q.p *= vol;
		Q.u *= Q.r;
		Q.v *= Q.r;
		Q.w *= Q.r;

		Q.add(D);
		Q.r = fmax(Q.r,smallr);
		Q.p = fmax(Q.p,smallp);

		G.C[ind].r = Q.r/vol;
		G.C[ind].p = Q.p/vol;
		G.C[ind].u = Q.u/Q.r;
		#if ndim>1
		G.C[ind].v = Q.v/Q.r;
		#endif
		#if ndim>2
		G.C[ind].w = Q.w/Q.r;
		#endif
	}

	return;
}
