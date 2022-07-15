__device__ double get_t_stop(double rho, double cs)
{
	return Stokes/(cs/sc_h)/rho;
}

__device__ Dust dust_flux(int i, int geom, double* xa, double* dx, double* dv, 
                          double* r, double* u, double* v, double* w, double dt, double us,
                          bool print=false)
{
	Dust F;
	State S;

	double r_par[4], u_par[4];

	double ul, ur;
	double xl, xr, tmp, q;
	double ql, qr;
	double r0, u0;
	double u_;

	/////////////////////////////////////////////////////////

	xl = xa[i];
	xr = xa[i+1];

	get_PRM_parameters(i, geom, xa, dx, dv, u, u_par);

	tmp = xl - fmin(u_par[1]-u_par[0],0.0)*dt;
	dimensionless_x(xl,tmp,xr,q,ql,qr);

	ur = get_PRM_aveL(geom, q, u_par, ql, qr);
	u_ = 0.5*dt*(ur-u[i])/(0.5*(xr-tmp));

	/////////////////////////////////////////////////////////

	tmp = xl - fmin(ur, 0.0)*dt;
	dimensionless_x(xl,tmp,xr,q,ql,qr);

	get_CON_parameters(i, geom, xa, dx, dv, r, r_par);

	r0 = get_CON_aveL(geom, q, r_par, ql, qr);
	u0 = get_PRM_aveL(geom, q, u_par, ql, qr);
	
	S.rr = r0*exp_lim(u_);
	S.ur = u0 + us;

	get_PRM_parameters(i, geom, xa, dx, dv, v, r_par);
	get_PRM_parameters(i, geom, xa, dx, dv, w, u_par);

	S.vr = get_PRM_aveL(geom, q, r_par, ql, qr);
	S.wr = get_PRM_aveL(geom, q, u_par, ql, qr);

	if (print) printf(" ur=%e\n r0=%e, u0=%e\n u_=%e\n",ur,r0,u0,u_);


	/////////////////////////////////////////////////////////

	xl = xa[i-1];
	xr = xa[i];

	get_PRM_parameters(i-1, geom, xa, dx, dv, u, u_par);

	tmp = xr - fmax(u_par[1]+u_par[2],0.0)*dt;
	dimensionless_x(xl,tmp,xr,q,ql,qr);

	ul = get_PRM_aveR(geom, q, u_par, ql, qr);
	u_ = 0.5*dt*(u[i-1]-ul)/(0.5*(tmp-xl));

	/////////////////////////////////////////////////////////

	tmp = xr - fmax(ul,0.0)*dt;
	dimensionless_x(xl,tmp,xr,q,ql,qr);

	get_CON_parameters(i-1, geom, xa, dx, dv, r, r_par);

	r0 = get_CON_aveR(geom, q, r_par, ql, qr);
	u0 = get_PRM_aveR(geom, q, u_par, ql, qr);

	S.rl = r0*exp_lim(u_);
	S.ul = u0 + us;

	get_PRM_parameters(i-1, geom, xa, dx, dv, v, r_par);
	get_PRM_parameters(i-1, geom, xa, dx, dv, w, u_par);

	S.vl = get_PRM_aveR(geom, q, r_par, ql, qr);
	S.wl = get_PRM_aveR(geom, q, u_par, ql, qr);

	if (print) printf(" ul=%e\n r0=%e, u0=%e\n u_=%e\n",ul,r0,u0,u_);

	/////////////////////////////////////////////////////////

	S.rl = fmin(r[i-1]*dt/(xa[i]-xa[i-1]),fmax(S.ul*S.rl,0.0));
	S.rr = fmax(  r[i]*dt/(xa[i]-xa[i+1]),fmin(S.ur*S.rr,0.0));

	F.r = S.rl + S.rr;
	F.u = S.rl*S.ul + S.rr*S.ur;
	F.v = S.rl*S.vl + S.rr*S.vr;
	F.w = S.rl*S.wl + S.rr*S.wr;

	return F;
}

__device__ Dust dust_flow(int geom, double* xa, double* dx, double* dv, double rad, double* r, double* u, double* v, double* w, double force, double dt)
{
	Dust Del;

	int imax = blockDim.x;
	int i = threadIdx.x;

	double us;

	#if recon_flag==2
	extern __shared__ double share[];
	double* flat = &share[0];
	int is = i + imax*threadIdx.y;
	flat[is] = 0.0;
	__syncthreads();
	#endif

	if (i>=npad && i<imax+1-npad)
	{
		us = 0.5*dt*force;
		if (geom>2) dt /= rad;

		Del = dust_flux(i, geom, xa, dx, dv, r, u, v, w, dt, us);
	}
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

__device__ Dust update_Dust(double x, double y, double z, Dust Q, Dust D)
{
	double r_min = fmax(Q.r*1.0e-10,D_G_ratio*smallr);
	double u_old = Q.u;
	double v_old = Q.v;
	double w_old = Q.w;

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
		Q.u = u_old;
		Q.v = v_old;
		Q.w = w_old;
		//printf("Error: negative density at %f %f %f\n",G.get_xc(i),G.get_yc(j),G.get_zc(k));
	}

	return Q;
}

__global__ void sweepx_dust_inplace(Grid G, Dust* in, Dust* out, double dt)
{
	__shared__ double xa[x_xthd+1], dx[x_xthd], xv[x_xthd];
	__shared__ double r[x_xthd*x_ydiv], u[x_xthd*x_ydiv], v[x_xthd*x_ydiv], w[x_xthd*x_ydiv];

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

	r[i+x_xthd*j] = in[ind].r;
	u[i+x_xthd*j] = in[ind].u;
	v[i+x_xthd*j] = in[ind].v;
	w[i+x_xthd*j] = in[ind].w;

	#if geomx == 1
	v[i+x_xthd*j] *= rad;
	#elif geomx == 2
	v[i+x_xthd*j] *= rad_cyl;
	w[i+x_xthd*j] *= rad;
	#endif
	__syncthreads();

	/////////////////////////////////////////////////////
	Dust Del;
	Del = dust_flow(geomx, xa, dx, xv, rad, &r[x_xthd*j], &u[x_xthd*j], &v[x_xthd*j], &w[x_xthd*j], 0.0, dt);
	Del.multiply(dt/xv[i]);

	if (i>=xpad && i<x_xthd-xpad)
	{
		#if geomx == 1
		Del.v /= rad;
		#elif geomx == 2
		Del.v /= rad_cyl;
		Del.w /= rad;
		#endif
		out[ind] = update_Dust(G.get_xc(idx), G.get_yc(idy), G.get_zc(idz), in[ind], Del);
	}
	return;
}

void sweepx_dust_inplace(Grid* dev, double dt)
{
	int nx, ny, nz;

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		sweepx_dust_inplace<<< dim3(nx/x_xdiv,ny/x_ydiv,nz/x_zdiv), dim3(x_xthd,x_ydiv,x_zdiv), 2*sizeof(double)*x_xthd*x_ydiv*x_zdiv, dev[n].stream >>>
		                   (dev[n], dev[n].CD, dev[n].TD, dt);
		dev[n].CT_D_change();
	}
	return;
}


__global__ void sweepy_dust_inplace(Grid G, Dust* in, Dust* out, double dt)
{
	__shared__ double ya[y_ythd+1], dy[y_ythd], yv[y_ythd];
	__shared__ double r[y_ythd*y_xdiv], u[y_ythd*y_xdiv], v[y_ythd*y_xdiv], w[y_ythd*y_xdiv];

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

	r[i+y_ythd*j] = in[ind].r;
	u[i+y_ythd*j] = in[ind].v - G.get_rot(idx,idz);
	v[i+y_ythd*j] = in[ind].w;
	w[i+y_ythd*j] = in[ind].u;

	#if geomy == 3 || geomy == 4
	u[i+y_ythd*j] -= rad_cyl*frame_omega;
	#endif
	__syncthreads();

	/////////////////////////////////////////////////////

	Dust Del;
	Del = dust_flow(geomy, ya, dy, yv, rad_cyl, &r[y_ythd*j], &u[y_ythd*j], &v[y_ythd*j], &w[y_ythd*j], 0.0, dt);
	Del.multiply(dt/yv[i]);

	#if geomy > 2
	Del.multiply(1.0/rad_cyl);
	#endif

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

		out[ind] = update_Dust(G.get_xc(idx), G.get_yc(idy), G.get_zc(idz), in[ind], Del);
	}

	return;
}

void sweepy_dust_inplace(Grid* dev, double dt)
{
	int nx, ny, nz;

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		sweepy_dust_inplace<<< dim3(ny/y_ydiv,nx/y_xdiv,nz/y_zdiv), dim3(y_ythd,y_xdiv,y_zdiv), 2*sizeof(double)*y_ythd*y_xdiv*y_zdiv, dev[n].stream >>>
		                   (dev[n], dev[n].CD, dev[n].TD, dt);
		dev[n].CT_D_change();
	}
	return;
}

__global__ void sweepz_dust_inplace(Grid G, Dust* in, Dust* out, double dt)
{
	__shared__ double za[z_zthd+1], dz[z_zthd], zv[z_zthd];
	__shared__ double r[z_zthd*z_xdiv], u[z_zthd*z_xdiv], v[z_zthd*z_xdiv], w[z_zthd*z_xdiv];

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

	r[i+z_zthd*j] = in[ind].r;
	u[i+z_zthd*j] = in[ind].w;
	v[i+z_zthd*j] = in[ind].u;
	w[i+z_zthd*j] = in[ind].v;

	#if geomz == 5
	w[i+z_zthd*j] *= rad_cyl;
	#endif
	__syncthreads();

	/////////////////////////////////////////////////////

	Dust Del;
	Del = dust_flow(geomz, za, dz, zv, rad, &r[z_zthd*j], &u[z_zthd*j], &v[z_zthd*j], &w[z_zthd*j], 0.0, dt);

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

		out[ind] = update_Dust(G.get_xc(idx), G.get_yc(idy), G.get_zc(idz), in[ind], Del);
	}

	return;
}

void sweepz_dust_inplace(Grid* dev, double dt)
{
	int nx, ny, nz;

	for (int n=0; n<ndev; n++)
	{
		cudaSetDevice(n);

		nx = dev[n].xres;
		ny = dev[n].yres;
		nz = dev[n].zres;

		sweepz_dust_inplace<<< dim3(nz/z_zdiv,nx/z_xdiv,ny/z_ydiv), dim3(z_zthd,z_xdiv,z_ydiv), 2*sizeof(double)*z_zthd*z_xdiv*z_ydiv, dev[n].stream >>>
		                   (dev[n], dev[n].CD, dev[n].TD, dt);
		dev[n].CT_D_change();
	}
	return;
}
