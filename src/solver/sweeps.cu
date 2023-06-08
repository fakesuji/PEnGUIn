__global__ void sweepx(Grid G, Cell* C, double dt)
{
	__shared__ double xa[x_xthd+1], dx[x_xthd], xv[x_xthd];
	__shared__ double r[x_xthd*x_ydiv], p[x_xthd*x_ydiv], u[x_xthd*x_ydiv], v[x_xthd*x_ydiv], w[x_xthd*x_ydiv];
	double force;

	int i = threadIdx.x;
	int idx = i + blockIdx.x*x_xdiv;

	int j = threadIdx.y;
	int idy = j + blockIdx.y*x_ydiv + ypad;

	int k = threadIdx.z;
	int idz = k + blockIdx.z*x_zdiv + zpad;

	int ind = G.get_ind(idx,idy,idz);

	if (j==0)
	{
		xa[i] = G.get_xa(idx);
		if (i==blockDim.x-1) xa[i+1] = G.get_xa(idx+1);
		xv[i] = G.get_xv(idx);
	}
	__syncthreads();

	r[i+x_xthd*j] = C[ind].r;
	p[i+x_xthd*j] = C[ind].p;
	u[i+x_xthd*j] = C[ind].u;
	v[i+x_xthd*j] = C[ind].v;
	w[i+x_xthd*j] = C[ind].w;
	if (i>0) force = 0.5*(G.fx[ind] + G.fx[G.get_ind(idx-1,idy,idz)]);

	__syncthreads();

	if (j==0) dx[i] = xa[i+1] - xa[i];

	double rad = 0.5*(xa[i+1]+xa[i]);
	#if geomx == 1
	v[i+x_xthd*j] *= rad;
	#elif geomx == 2
	double rad_cyl = rad * sin(G.get_zc(idz));
	v[i+x_xthd*j] *= rad_cyl;
	w[i+x_xthd*j] *= rad;
	#endif
	__syncthreads();

//if (idy==2) printf("%f %f %f\n",rad,r[y_ythd*j+i],p[y_ythd*j+i]);

	/////////////////////////////////////////////////////
	Cell Del;
	Del =   riemann(geomx, xa, dx, xv, rad, &r[x_xthd*j], &p[x_xthd*j], &u[x_xthd*j], &v[x_xthd*j], &w[x_xthd*j], force, dt);
//if (idy==ypad && i>=xpad && i<x_xthd-xpad && idx>270 && idx<290) printf("%f %f: %f, %f, %e, %f, %f, %e\n",xa[i],xv[i+1]-xv[i],r[i+x_xthd*j],p[i+x_xthd*j],u[i+x_xthd*j],force,G.get_rot(idx,idz),Del.r);

	Del.multiply(G.get_yv(idy)*G.get_zv(idz));

	if (i>=xpad && i<x_xthd-xpad)
	{
		#if geomx == 1
		Del.v /= rad;
		#elif geomx == 2
		Del.v /= rad_cyl;
		Del.w /= rad;
		#endif
		G.T[ind].add(Del);
	}
	//if (idx==12 && G.get_j_shf(idx,idy,idz)==17) printf("sweepx: rad=%f azi=%f\n",G.get_xc(idx),G.get_yc(idy));

	return;
}


__global__ void sweepy(Grid G, Cell* C, double dt)
{
	__shared__ double ya[y_ythd+1], dy[y_ythd], yv[y_ythd];
	__shared__ double r[y_ythd*y_xdiv], p[y_ythd*y_xdiv], u[y_ythd*y_xdiv], v[y_ythd*y_xdiv], w[y_ythd*y_xdiv];
	double force;

	int i = threadIdx.x;
	int idy = i + blockIdx.x*y_ydiv;

	int j = threadIdx.y;
	int idx = j + blockIdx.y*y_xdiv + xpad;

	int k = threadIdx.z;
	int idz = k + blockIdx.z*x_zdiv + zpad;

	int ind = idx + G.xarr*idy + G.xarr*G.yarr*idz;

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

	r[i+y_ythd*j] = C[ind].r;
	p[i+y_ythd*j] = C[ind].p;
	u[i+y_ythd*j] = C[ind].v - G.get_rot(idx,idz);
	v[i+y_ythd*j] = C[ind].w;
	w[i+y_ythd*j] = C[ind].u;

	#if geomy == 3 || geomy == 4
	u[i+y_ythd*j] -= rad*frame_omega;
	#endif

	if (i>0) force = 0.5*(G.fy[ind] + G.fy[G.get_ind(idx,idy-1,idz)]);
	__syncthreads();

	/////////////////////////////////////////////////////

	Cell Del;
	Del = riemann(geomy, ya, dy, yv, rad, &r[y_ythd*j], &p[y_ythd*j], &u[y_ythd*j], &v[y_ythd*j], &w[y_ythd*j], force, dt);

	Del.multiply(G.get_xv(idx)*G.get_zv(idz));

	#if geomy > 2
	Del.multiply(1.0/rad);
	#endif

	if (i>=ypad && i<y_ythd-ypad)
	{
		G.T[ind].r += Del.r;
		G.T[ind].p += Del.p;
		G.T[ind].u += Del.w;
		#if geomy == 3 || geomy == 4
		G.T[ind].v += Del.u + (G.get_rot(idx,idz) + rad*frame_omega)*Del.r;
		#else
		G.T[ind].v += Del.u + G.get_rot(idx,idz)*Del.r;
		#endif
		G.T[ind].w += Del.v;
	}

	return;
}

__global__ void sweepz(Grid G, Cell* C, double dt)
{
	__shared__ double za[z_zthd+1], dz[z_zthd], zv[z_zthd];
	__shared__ double r[z_zthd*z_xdiv], p[z_zthd*z_xdiv], u[z_zthd*z_xdiv], v[z_zthd*z_xdiv], w[z_zthd*z_xdiv];
	double force;

	int i = threadIdx.x;
	int idz = i + blockIdx.x*z_zdiv;

	int j = threadIdx.y;
	int idx = j + blockIdx.y*z_xdiv + xpad;

	int k = threadIdx.z;
	int idy = k + blockIdx.z*z_ydiv + ypad;

	int ind = G.get_ind(idx,idy,idz);

	if (j==0)
	{
		za[i] = G.get_za(idz);
		if (i==blockDim.x-1) za[i+1] = G.get_za(idz+1);
		zv[i] = G.get_zv(idz);
	}
	__syncthreads();

	if (j==0) dz[i] = za[i+1] - za[i];
	
	r[i+z_zthd*j] = C[ind].r;
	p[i+z_zthd*j] = C[ind].p;
	u[i+z_zthd*j] = C[ind].w;
	v[i+z_zthd*j] = C[ind].u;
	w[i+z_zthd*j] = C[ind].v;
	if (i>0) force = 0.5*(G.fz[ind] + G.fz[G.get_ind(idx,idy,idz-1)]);
	__syncthreads();

	double rad;
	rad = G.get_xc(idx);
	#if geomz == 5
	double rad_cyl;
	rad_cyl = rad * sin(0.5*(za[i+1]+za[i]));
	w[i+z_zthd*j] *= rad_cyl;
	#endif
	__syncthreads();

	/////////////////////////////////////////////////////

	Cell Del;
	Del = riemann(geomz, za, dz, zv, rad, &r[z_zthd*j], &p[z_zthd*j], &u[z_zthd*j], &v[z_zthd*j], &w[z_zthd*j], force, dt);
	Del.multiply(G.get_xv(idx)*G.get_yv(idy));

	#if geomz > 2
	Del.multiply(1.0/rad);
	#endif

	if (i>=zpad && i<z_zthd-zpad)
	{
		#if geomz == 5
		Del.w /= rad_cyl;
		#endif
		G.T[ind].r += Del.r;
		G.T[ind].p += Del.p;
		G.T[ind].u += Del.v;
		G.T[ind].v += Del.w;
		G.T[ind].w += Del.u;
	}

	return;
}

__device__ Cell update_Cell(double x, double y, double z, Cell Q, Cell D)
{
	double r_min = fmax(Q.r*1.0e-10,smallr);
	double p_min = fmax(Q.p*1.0e-10,smallp);
	double u_old = Q.u;
	double v_old = Q.v;
	double w_old = Q.w;

	Q.p = get_energy(Q.r,Q.p,Q.u,Q.v,Q.w);
	Q.p *= Q.r;
	Q.u *= Q.r;
	Q.v *= Q.r;
	Q.w *= Q.r;

	Q.add(D);
	
	Q.u /= Q.r;
	Q.v /= Q.r;
	Q.w /= Q.r;

	if (Q.r<r_min)
	{
		Q.r = r_min;
		Q.p = p_min;
		Q.u = u_old;
		Q.v = v_old;
		Q.w = w_old;
		//printf("Error: negative density at %f %f %f\n",G.get_xc(i),G.get_yc(j),G.get_zc(k));
	}
	else
	{
		#if EOS_flag == 2
		#if internal_e_flag==0
		Q.p = fmax(Q.p*gamm - gamm*Q.r*(Q.u*Q.u+Q.v*Q.v+Q.w*Q.w)/2.0,p_min);
		#else
		Q.p = fmax(Q.p*gamm,p_min);
		#endif
		#elif EOS_flag == 0
		Q.p = get_cs2(x, y, z)*Q.r;
		#endif
	}

	return Q;
}

__global__ void sweepx_inplace(Grid G, Cell* C, Cell* out, double dt)
{
	__shared__ double xa[x_xthd+1], dx[x_xthd], xv[x_xthd];
	__shared__ double r[x_xthd*x_ydiv], p[x_xthd*x_ydiv], u[x_xthd*x_ydiv], v[x_xthd*x_ydiv], w[x_xthd*x_ydiv];

	int i = threadIdx.x;
	int idx = i + blockIdx.x*x_xdiv;

	int j = threadIdx.y;
	int idy = j + blockIdx.y*x_ydiv + ypad;

	int k = threadIdx.z;
	int idz = k + blockIdx.z*x_zdiv + zpad;

	int ind = G.get_ind(idx,idy,idz);

	if (j==0)
	{
		xa[i] = G.get_xa(idx);
		if (i==blockDim.x-1) xa[i+1] = G.get_xa(idx+1);
		xv[i] = G.get_xv(idx);
	}
	__syncthreads();

	if (j==0) dx[i] = xa[i+1] - xa[i];

	double rad = 0.5*(xa[i+1]+xa[i]);
	#if geomx == 2
	double rad_cyl = rad * sin(G.get_zc(idz));
	#endif

	r[i+x_xthd*j] = C[ind].r;
	p[i+x_xthd*j] = C[ind].p;
	u[i+x_xthd*j] = C[ind].u;
	v[i+x_xthd*j] = C[ind].v;
	w[i+x_xthd*j] = C[ind].w;

	#if geomx == 1
	v[i+x_xthd*j] *= rad;
	#elif geomx == 2
	v[i+x_xthd*j] *= rad_cyl;
	w[i+x_xthd*j] *= rad;
	#endif
	__syncthreads();

	/////////////////////////////////////////////////////
	Cell Del;
	Del = riemann(geomx, xa, dx, xv, 1.0, &r[x_xthd*j], &p[x_xthd*j], &u[x_xthd*j], &v[x_xthd*j], &w[x_xthd*j], 0.0, dt);
	Del.multiply(dt/xv[i]);

	if (i>=xpad && i<x_xthd-xpad)
	{
		#if geomx == 1
		Del.v /= rad;
		#elif geomx == 2
		Del.v /= rad_cyl;
		Del.w /= rad;
		#endif
<<<<<<< HEAD
/*
		if (idx<5 && idy==100) printf("R=%f; Del_r=%e\n",rad,Del.r);
		if (idx==10 && idy==100) 
=======

		//if (idx<5 && idy==100) printf("R=%f; Del_r=%e\n",rad,Del.r);
		/*if (idx==2)// && idy==100) 
>>>>>>> methods
		{
			printf("%f, %f, %f, %f, %f\n",xa[i-2],xa[i-1],xa[i],xa[i+1],xa[i+2]);
			printf("%e, %e, %e, %e, %e\n",r[i+x_xthd*j-2],r[i+x_xthd*j-1],r[i+x_xthd*j],r[i+x_xthd*j+1],r[i+x_xthd*j+2]);
			printf("%e, %e, %e, %e, %e\n",u[i+x_xthd*j-2],u[i+x_xthd*j-1],u[i+x_xthd*j],u[i+x_xthd*j+1],u[i+x_xthd*j+2]);
			printf("%e, %e, %e\n\n",Del.r,Del.p,Del.u);
		}
		*/

		out[ind] = update_Cell(G.get_xc(idx), G.get_yc(idy), G.get_zc(idz), C[ind], Del);
	}
	return;
}

void sweepx_inplace(Grid* dev, double dt)
{
	int nx, ny, nz;
	Cell* Cp;

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		sweepx_inplace<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 2*sizeof(double)*x_xthd*x_ydiv*x_zdiv, dev[n].stream >>>
		      (dev[n], dev[n].C, dev[n].T, dt);

		Cp = dev[n].C;
		dev[n].C = dev[n].T;
		dev[n].T = Cp;
	}
	return;
}


__global__ void sweepy_inplace(Grid G, Cell* C, Cell* out, double dt)
{
	__shared__ double ya[y_ythd+1], dy[y_ythd], yv[y_ythd];
	__shared__ double r[y_ythd*y_xdiv], p[y_ythd*y_xdiv], u[y_ythd*y_xdiv], v[y_ythd*y_xdiv], w[y_ythd*y_xdiv];

	int i = threadIdx.x;
	int idy = i + blockIdx.x*y_ydiv;

	int j = threadIdx.y;
	int idx = j + blockIdx.y*y_xdiv + xpad;

	int k = threadIdx.z;
	int idz = k + blockIdx.z*x_zdiv + zpad;

	int ind = idx + G.xarr*idy + G.xarr*G.yarr*idz;

	if (j==0)
	{
		ya[i] = G.get_ya(idy);
		if (i==blockDim.x-1) ya[i+1] = G.get_ya(idy+1);
		yv[i] = G.get_yv(idy);
	}
	__syncthreads();

	if (j==0) dy[i] = ya[i+1] - ya[i];

	double rad_cyl;	
	#if geomy == 3
	rad_cyl = G.get_xc(idx);
	#elif geomy == 4
	rad_cyl = G.get_xc(idx)*sin(G.get_zc(idz));
	#else
	rad_cyl = 1.0;
	#endif

	r[i+y_ythd*j] = C[ind].r;
	p[i+y_ythd*j] = C[ind].p;
	u[i+y_ythd*j] = C[ind].v - G.get_rot(idx,idz);
	v[i+y_ythd*j] = C[ind].w;
	w[i+y_ythd*j] = C[ind].u;

	#if geomy == 3 || geomy == 4
	u[i+y_ythd*j] -= rad_cyl*frame_omega;
	#endif
	__syncthreads();

	/////////////////////////////////////////////////////

	Cell Del;
	Del = riemann(geomy, ya, dy, yv, rad_cyl, &r[y_ythd*j], &p[y_ythd*j], &u[y_ythd*j], &v[y_ythd*j], &w[y_ythd*j], 0.0, dt);

	Del.multiply(dt/yv[i]/rad_cyl);

	double tmp;

	if (i>=ypad && i<y_ythd-ypad)
	{
		#if geomy == 3 || geomy == 4
		Del.u += (G.get_rot(idx,idz) + rad_cyl*frame_omega)*Del.r;
		#else
		Del.u += G.get_rot(idx,idz)*Del.r;
		#endif

		tmp = Del.v;
		Del.v = Del.u;
		Del.u = Del.w;
		Del.w = tmp;

		out[ind] = update_Cell(G.get_xc(idx), G.get_yc(idy), G.get_zc(idz), C[ind], Del);
	}

	return;
}

void sweepy_inplace(Grid* dev, double dt)
{
	int nx, ny, nz;
	Cell* Cp;

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		sweepy_inplace<<< dim3(ny/y_ydiv,nx/y_xdiv,nz/y_zdiv), dim3(y_ythd,y_xdiv,y_zdiv), 2*sizeof(double)*y_ythd*y_xdiv*y_zdiv, dev[n].stream >>>
		      (dev[n], dev[n].C, dev[n].T, dt);

		Cp = dev[n].C;
		dev[n].C = dev[n].T;
		dev[n].T = Cp;
	}
	return;
}

__global__ void sweepz_inplace(Grid G, Cell* C, Cell* out, double dt)
{
	__shared__ double za[z_zthd+1], dz[z_zthd], zv[z_zthd];
	__shared__ double r[z_zthd*z_xdiv], p[z_zthd*z_xdiv], u[z_zthd*z_xdiv], v[z_zthd*z_xdiv], w[z_zthd*z_xdiv];

	int i = threadIdx.x;
	int idz = i + blockIdx.x*z_zdiv;

	int j = threadIdx.y;
	int idx = j + blockIdx.y*z_xdiv + xpad;

	int k = threadIdx.z;
	int idy = k + blockIdx.z*z_ydiv + ypad;

	int ind = G.get_ind(idx,idy,idz);

	if (j==0)
	{
		za[i] = G.get_za(idz);
		if (i==blockDim.x-1) za[i+1] = G.get_za(idz+1);
		zv[i] = G.get_zv(idz);
	}
	__syncthreads();

	if (j==0) dz[i] = za[i+1] - za[i];
	
	#if geomz == 5
	double rad = G.get_xc(idx);
	double rad_cyl = rad * sin(0.5*(za[i+1]+za[i]));
	#else
	double rad = 1.0;
	#endif

	r[i+z_zthd*j] = C[ind].r;
	p[i+z_zthd*j] = C[ind].p;
	u[i+z_zthd*j] = C[ind].w;
	v[i+z_zthd*j] = C[ind].u;
	w[i+z_zthd*j] = C[ind].v;

	#if geomz == 5
	w[i+z_zthd*j] *= rad_cyl;
	#endif
	__syncthreads();

	/////////////////////////////////////////////////////

	Cell Del;
	Del = riemann(geomz, za, dz, zv, rad, &r[z_zthd*j], &p[z_zthd*j], &u[z_zthd*j], &v[z_zthd*j], &w[z_zthd*j], 0.0, dt);

	Del.multiply(dt/zv[i]);

	#if geomz == 5
	Del.multiply(1.0/rad);
	#endif

	double tmp;

	if (i>=zpad && i<z_zthd-zpad)
	{
		#if geomz == 5
		Del.w /= rad_cyl;
		#endif

		tmp = Del.w;
		Del.w = Del.u;
		Del.u = Del.v;
		Del.v = tmp;

		out[ind] = update_Cell(G.get_xc(idx), G.get_yc(idy), G.get_zc(idz), C[ind], Del);
	}

	return;
}

void sweepz_inplace(Grid* dev, double dt)
{
	int nx, ny, nz;
	Cell* Cp;

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		sweepz_inplace<<< dim3(nz/z_zdiv,nx/z_xdiv,ny/z_ydiv), dim3(z_zthd,z_xdiv,z_ydiv), 2*sizeof(double)*z_zthd*z_xdiv*z_ydiv, dev[n].stream >>>
		      (dev[n], dev[n].C, dev[n].T, dt);

		Cp = dev[n].C;
		dev[n].C = dev[n].T;
		dev[n].T = Cp;
	}
	return;
}

__device__ void get_g_forces(Grid G, Cell T, Dust D, double ts, double xc, double yc, double zc, int ind, double &fx, double &fy, double &fz, Cell C)
{
	#if ndim > 2
	fz = get_fz(xc,yc,zc,T.u,T.v,T.w);
	#ifdef visc_flag
	fz += G.vis_tensor[ndim*ind+2]/T.r;
	#endif
	#endif

	#if ndim > 1
	fy = get_fy(xc,yc,zc,T.u,T.v,T.w);
	#ifdef visc_flag
	fy += G.vis_tensor[ndim*ind+1]/T.r;
	#endif
	#endif

	fx = get_fx(xc,yc,zc,T.u,T.v,T.w);
	#ifdef visc_flag
	fx += G.vis_tensor[ndim*ind+0]/T.r;
	#endif
	return;
}

__device__ void get_d_forces(Grid G, Cell T, Dust D, double ts, double xc, double yc, double zc, int ind, double &fx, double &fy, double &fz, Dust Q)
{
	#if ndim > 2
	fz = get_fz(xc,yc,zc,D.u,D.v,D.w);
	fz = (fz*ts + T.w - Q.w);
	#endif

	#if ndim > 1
	fy = get_fy(xc,yc,zc,D.u,D.v,D.w);
	fy = (fy*ts + T.v - Q.v);
	#endif

	fx = get_fx(xc,yc,zc,D.u,D.v,D.w);
	fx = (fx*ts + T.u - Q.u);

	return;
}


__global__ void apply_source_terms(Grid G, Cell* in, Dust* in_d, Cell* out, Dust* out_d, double x_dt, double mdt, double dt)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.y + blockIdx.y*blockDim.y;
	int k = threadIdx.z + blockIdx.z*blockDim.z;

	Cell C, T;
	Dust D;
	double xc,yc,zc;
	double fx = 0.0;
	double fy = 0.0;
	double fz = 0.0;
	double gx = 0.0;
	double gy = 0.0;
	double gz = 0.0;
	double ts = 0.0;
	int ind;

	#ifdef dust_flag
	double fx_d = 0.0;
	double fy_d = 0.0;
	double fz_d = 0.0;
	double dtau;
	Dust Q;
	#endif

	if (i>=xpad-1 && i<G.xarr-xpad+1)
	if (j>=ypad-1 && j<G.yarr-ypad+1)
	if (k>=zpad-1 && k<G.zarr-zpad+1)
	{		
		ind = G.get_ind(i,j,k);
		C.copy(in[ind]);
		T.copy(C);
		#ifdef dust_flag
		Q.copy(in_d[ind]);
		D.copy(Q);
		#endif

		////////////////////////////////////////////////

		xc = G.get_xc(i);
		#if ndim > 1
		yc = G.get_yc(j) + G.get_rot(i,k)*x_dt;
		#else
		yc = 0.0;
		#endif
		#if ndim > 2
		zc = G.get_zc(k);
		#else
		zc = 0.0;
		#endif

		gx = get_gx(xc, yc, zc, G.planets);
		gy = get_gy(xc, yc, zc, G.planets);
		gz = get_gz(xc, yc, zc, G.planets);

		////////////////////////////////////////////////

		#ifdef dust_flag
		ts = get_t_stop(T.r,sqrt(gam*T.p/T.r));
		get_d_forces(G, T, D, ts, xc, yc, zc, ind, fx_d, fy_d, fz_d, Q);
		#endif
		get_g_forces(G, T, D, ts, xc, yc, zc, ind, fx, fy, fz, C);

		////////////////////////////////////////////////

		if (mdt!=0.0)
		{
			T.u = C.u + fx*mdt + gx*mdt;
			T.v = C.v + fy*mdt + gy*mdt;
			T.w = C.w + fz*mdt + gz*mdt;
			#ifdef cool_flag
			T.p += cooling(xc,yc,zc,T.p,T.r,mdt);
			#endif

			#ifdef dust_flag
			dtau = 1.0-exp(-mdt/ts);
	
			D.u = Q.u + fx_d*dtau + gx*ts*dtau;
			D.v = Q.v + fy_d*dtau + gy*ts*dtau;
			D.w = Q.w + fz_d*dtau + gz*ts*dtau;
			#endif

			////////////////////////////////////////////////
	
			#ifdef dust_flag
			ts = get_t_stop(T.r,sqrt(gam*T.p/T.r));
			get_d_forces(G, T, D, ts, xc, yc, zc, ind, fx_d, fy_d, fz_d, Q);
			#endif
			get_g_forces(G, T, D, ts, xc, yc, zc, ind, fx, fy, fz, C);
		}

		////////////////////////////////////////////////

		out[ind].u = C.u + fx*dt + gx*dt;
		out[ind].v = C.v + fy*dt + gy*dt;
		out[ind].w = C.w + fz*dt + gz*dt;
		#ifdef cool_flag
		out[ind].p += cooling(xc,yc,zc,out[ind].p,out[ind].r,dt);
		#endif

		#ifdef dust_flag
		dtau = 1.0-exp(-dt/ts);
		out_d[ind].u = Q.u + fx_d*dtau + gx*ts*dtau;
		out_d[ind].v = Q.v + fy_d*dtau + gy*ts*dtau;
		out_d[ind].w = Q.w + fz_d*dtau + gz*ts*dtau;
		#endif
	}

	return;
}

void apply_source_terms(Grid* dev, double x_dt, double mdt, double dt)
{
	int mx, my, mz;
	int bsz = 32;

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		mx = dev[n].xarr;
		my = dev[n].yarr;
		mz = dev[n].zarr;

		apply_source_terms<<< dim3((mx+bsz-1)/bsz,my,mz), bsz, 0, dev[n].stream >>> (dev[n], dev[n].C, dev[n].CD, dev[n].T, dev[n].TD, x_dt, mdt, dt);
	}
	return;
}

void apply_source_terms_inplace(Grid* dev, double x_dt, double mdt, double dt)
{
	int mx, my, mz;
	int bsz = 32;

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		mx = dev[n].xarr;
		my = dev[n].yarr;
		mz = dev[n].zarr;

		apply_source_terms<<< dim3((mx+bsz-1)/bsz,my,mz), bsz, 0, dev[n].stream >>> (dev[n], dev[n].C, dev[n].CD, dev[n].C, dev[n].CD, x_dt, mdt, dt);
	}
	return;
}
