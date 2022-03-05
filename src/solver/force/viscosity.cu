__device__ void derivs(double &first, double &second, double x2, double x3, double u2, double u3)
{
  double h2, h3, A, B;

  h2 = u2/x2;
  h3 = u3/x3;

  A = (   h3-   h2)/(x3-x2);
  B = (x3*h2-x2*h3)/(x3-x2);

  second = 2.0*A;
  first  = second*x2 + B;
  
  return;
}

__device__ void derivs(double &first, double &second, int &n, double* xf, double* u)
{
  double u1 = u[n-1];
  double x1 = xf[n-1];
  derivs(first, second, xf[n]-x1, xf[n+1]-x1, u[n]-u1, u[n+1]-u1);

  return;
}

__device__ void derivs(double &first, double &second, int &n, double &x2, double &x3, double* u)
{
  double u1 = u[n-1];
  derivs(first, second, x2, x3, u[n]-u1, u[n+1]-u1);

  return;
}

__device__ double deriv_1st(double u2, double u3, double x2, double x3)
{
  double h2, h3;

  h2 = u2/x2;
  h3 = u3/x3;
  
  return (2.0*x2*(h3-h2) + (x3*h2-x2*h3))/(x3-x2);
}

__device__ double deriv_2nd(double u2, double u3, double x2, double x3)
{
  double h2, h3;

  h2 = u2/x2;
  h3 = u3/x3;
  
  return 2.0*(h3-h2)/(x3-x2);
}

__device__ double cross_deriv(double x1y1, double x2y1, double x1y2, double x2y2, double dx, double dy)
{
  return (x2y2 + x1y1 + x1y2 + x2y1)/(dx*dy);
}

//================================================================================
__device__ void vis_xx_bound(double &fx, double &fy, int i, double *r, double *u, double *v, double *w, double *x)
{
  double rho, R, d1, d2, dn;
  double n1, n2, n3, u1, u2, u3, v1, v2, v3, x1, x2, x3;

  rho = 0.5*(r[i]+r[i-1]);
  R = x[i];

  x1 = x[i-1];
  x2 = R;
  x3 = x[i+1];

  u1 = 0.5*(u[i-1]+u[i-2]);
  u2 = 0.5*(u[i]+u[i-1]);
  u3 = 0.5*(u[i+1]+u[i]);

  v1 = 0.5*(v[i-1]+v[i-2])/x1/x1;
  v2 = 0.5*(v[i]+v[i-1])/R/R;
  v3 = 0.5*(v[i+1]+v[i])/x3/x3;

  n1 = get_nu(x1)*0.5*(r[i-1]+r[i-2]);
  n2 = get_nu(x2)*rho;
  n3 = get_nu(x3)*0.5*(r[i+1]+r[i]);

  x2 -= x1;
  x3 -= x1;

  u2 -= u1;
  u3 -= u1;

  v2 -= v1;
  v3 -= v1;

  n2 -= n1;
  n3 -= n1;

  dn = deriv_1st(n2, n3, x2, x3);
  n2 += n1;

  derivs(d1, d2, x2, x3, u2, u3);
  u2 += u1;
  fx = (2.0*n2*d2 + 2.0*dn*d1 + 2.0*n2*d1/R - 2.0*n2*u2/R/R)/rho;

  derivs(d1, d2, x2, x3, v2, v3);
  fy = R*(R*n2*d2 + R*dn*d1 + 3.0*n2*d1)/rho;

  return;
}

__device__ void vis_xx(double &fx, double &fy, int i, double *r, double *u, double *v, double *w, double *x)
{
  double rho, R, d1, d2, dn;
  double n1, n2, n3, u1, u2, u3, v1, v2, v3, x1, x2, x3;

  rho = r[i];
  R = x[i];

  x1 = x[i-1];
  x2 = R;
  x3 = x[i+1];

  u1 = u[i-1];
  u2 = u[i];
  u3 = u[i+1];  

  v1 = v[i-1]/x1/x1;
  v2 = v[i]/R/R;
  v3 = v[i+1]/x3/x3;

  n1 = get_nu(x1)*r[i-1];
  n2 = get_nu(x2)*rho;
  n3 = get_nu(x3)*r[i+1];

  x2 -= x1;
  x3 -= x1;

  u2 -= u1;
  u3 -= u1;

  v2 -= v1;
  v3 -= v1;

  n2 -= n1;
  n3 -= n1;

  dn = deriv_1st(n2, n3, x2, x3);
  n2 += n1;

  derivs(d1, d2, x2, x3, u2, u3);
  u2 += u1;
  fx = (2.0*n2*d2 + 2.0*dn*d1 + 2.0*n2*d1/R - 2.0*n2*u2/R/R)/rho;

  derivs(d1, d2, x2, x3, v2, v3);
  fy = R*(R*n2*d2 + R*dn*d1 + 3.0*n2*d1)/rho;

  return;
}

//================================================================================

__device__ void device_viscosity_r(int i, int lim, double dt, double time, double R, double *x, double *cell_r, double *u, double *v, double *w, double *n)
{
  double rho, nr, d1, d2, dn;
  double ddu, ddv;
  #if ndim == 3
  double ddw;
  #endif

  if (i>=0 && i<lim)
  {
    rho = cell_r[i];
    nr = get_nu(R, time)*rho;
    n[i] = nr;

    v[i] /= R*R;
  }
  __syncthreads();

  if (i>=n_pad && i<lim-n_pad)
  {
    derivs(dn, d2, i, x, n);

    derivs(d1, d2, i, x, u);
    ddu = 2.0*nr*d2 + 2.0*dn*d1 + 2.0*nr*d1/R - 2.0*nr*u[i]/R/R;

    derivs(d1, d2, i, x, v);
    ddv = R*nr*d2 + R*dn*d1 + 3.0*nr*d1;

    #if ndim == 3
    derivs(d1, d2, i, x, w);
    ddw = nr*d2 + dn*d1 + nr*d1/R;
    #endif

    ddu = dt*ddu/rho;
    ddv = dt*ddv/rho;
    #if ndim == 3
    ddw = dt*ddw/rho;
    #endif
  }
  __syncthreads();

  if (i>=0 && i<lim) v[i] *= R*R;

  if (i>=n_pad && i<lim-n_pad)
  {
    u[i] += ddu;
    v[i] += ddv*R;
    #if ndim == 3
    w[i] += ddw;
    #endif
  }
  __syncthreads();

  return;
}

__device__ void device_viscosity_p(int j, int lim, double dt, double time, double R, double *x, double *cell_r, double *u, double *v, double *w, double *n)
{
  double rho, nr, d1, d2, dn;
  double ddu, ddv;
  #if ndim == 3
  double ddw;
  #endif

  if (j>=0 && j<lim)
  {
    rho = cell_r[j];
    nr = get_nu(R, time)*rho;
    n[j] = nr;

    x[j] *= R;
  }
  __syncthreads();

  if (j>=n_pad && j<lim-n_pad)
  {
    derivs(dn, d2, j, x, n);

    derivs(d1, d2, j, x, v);
    ddv = 2.0*nr*d2 + 2.0*dn*d1;
    ddu =-2.0*nr*d1/R;

    #if ndim == 3
    derivs(d1, d2, j, x, w);
    ddw = nr*d2 + dn*d1;
    #endif

    derivs(d1, d2, j, x, u);
    ddu += nr*d2 + dn*d1;
    ddv += 2.0*nr*d1/R + 2.0*dn*u[j]/R;

    ddv = dt*ddv/rho;
    #if ndim == 3
    ddw = dt*ddw/rho;
    #endif
    ddu = dt*ddu/rho;
  }
  __syncthreads();

  if (j>=0 && j<lim) x[j] /= R;

  if (j>=n_pad && j<lim-n_pad)
  {
    v[j] += ddv;
    #if ndim == 3
    w[j] += ddw;
    #endif
    u[j] += ddu;
  }
  __syncthreads();

  return;
}

__device__ void device_viscosity_z(int k, int lim, double dt, double time, double R, double *x, double *cell_r, double *u, double *v, double *w, double *n)
{
  double rho, nr, d1, d2, dn;
  double ddu, ddv, ddw;

  if (k>=0 && k<lim)
  {
    rho = cell_r[k];
    nr = get_nu(R, time)*rho;
    n[k] = nr;

    v[k] /= R;
  }
  __syncthreads();

  if (k>0 && k<arrsize-1)
  {
    derivs(dn, d2, k, x, n);
    derivs(d1, d2, k, x, w);
    ddw = 2.0*nr*d2 + 2.0*dn*d1; 

    derivs(d1, d2, k, x, u);
    ddu = nr*d2 + dn*d1;
    ddw+= nr*d1/R;

    derivs(d1, d2, k, x, v);
    ddv = nr*d2 + dn*d1;

    ddw = dt*ddw/rho;
    ddu = dt*ddu/rho;
    ddv = dt*ddv/rho;
  }
  __syncthreads();

  if (k>=0 && k<lim) v[k] *= R;

  if (k>0 && k<arrsize-1)
  {
    w[k] += ddw;
    u[k] += ddu;
    v[k] += ddv;
  }
  __syncthreads();

  return;
}

__device__ void load_cross(int istart, int iblk, ring_geometry *ygrid, ring *rings, ring *lft, ring *rgh, ring *udr, ring *top,
                           double &rad, double &dr, double &dphi, double &dz, double &rho, double &u, double &v, double &w, 
                           double &ur, double &vr, double &wr, double &ul, double &vl, double &wl, double &uf, double &vf, double &wf,
                           double &ub, double &vb, double &wb, double &ut, double &vt, double &wt, double &ud, double &vd, double &wd)
{
  int n = threadIdx.x + blockDim.x*blockIdx.x;

  int i = n/jmax;
  int j = n - i*jmax;
  int k = blockIdx.y;
  i--;
  #if ndim>2
  k--;
  #endif

  ring *cur;

  if (i>=iblk)       cur = &rgh[i- iblk + n_pad*k];
  else if (i<0)      cur = &lft[i+n_pad + n_pad*k];
  else if (k>=kmax)  cur = &top[k- kmax + n_pad*(i+istart)];
  else if (k<0)      cur = &udr[k+n_pad + n_pad*(i+istart)];
  else               cur = &rings[i + iblk*k];

  int jj = j_bound(j - (*cur).rot_j);

  rho = (*cur).C[jj].r;
  u = (*cur).C[jj].u;
  v = (*cur).C[jj].v;
  w = (*cur).C[jj].w;
  rad = (*cur).xc;

  #if ndim>1

  j  = j_bound(j+1);
  jj = j_bound(j - (*cur).rot_j);

  ur = (*cur).C[jj].u;
  vr = (*cur).C[jj].v;
  wr = (*cur).C[jj].w;
  dphi = (*ygrid).yc[j];
  j  = j_bound(j-1);

  ///////////////////////////////////////

  j  = j_bound(j-1);
  jj = j_bound(j - (*cur).rot_j);

  ul = (*cur).C[jj].u;
  vl = (*cur).C[jj].v;
  wl = (*cur).C[jj].w;
  dphi-= (*ygrid).yc[j];
  j  = j_bound(j+1);

  #endif

  i++;
  if (k==kmax) 
  {
    if (i==iblk) cur = &top[0 + n_pad*(i+istart-1)];
    else         cur = &top[0 + n_pad*(i+istart)];
  }
  else if (k==-1)
  {
    if (i==iblk) cur = &udr[n_pad-1 + n_pad*(i+istart-1)];
    else         cur = &udr[n_pad-1 + n_pad*(i+istart)];
  }
  else
  {
    if (i>=iblk) cur = &rgh[i-iblk + n_pad*k];
    else         cur = &rings[i + iblk*k];
  }
  jj = j_bound(j - (*cur).rot_j);

  uf = (*cur).C[jj].u;
  vf = (*cur).C[jj].v;
  wf = (*cur).C[jj].w;
  dr = (*cur).xc;
  i--;

  ///////////////////////////////////////

  i--;
  if (k==kmax) 
  {
    if (i==-1) cur = &top[0 + n_pad*(i+istart+1)];
    else       cur = &top[0 + n_pad*(i+istart)];
  }
  else if (k==-1)
  {
    if (i==-1) cur = &udr[n_pad-1 + n_pad*(i+istart+1)];
    else       cur = &udr[n_pad-1 + n_pad*(i+istart)];
  }
  else
  {
    if (i<=-1) cur = &lft[i+n_pad + n_pad*k];
    else       cur = &rings[i + iblk*k];
  }
  jj = j_bound(j - (*cur).rot_j);

  ub = (*cur).C[jj].u;
  vb = (*cur).C[jj].v;
  wb = (*cur).C[jj].w;
  dr-= (*cur).xc;;
  i++;

  #if ndim>2

  k++;
  if (i==iblk) 
  {
    if (k==kmax) cur = &rgh[0 + n_pad*(k-1)];
    else         cur = &rgh[0 + n_pad*k];
  }
  else if (i==-1)
  {
    if (k==kmax) cur = &lft[n_pad-1 + n_pad*(k-1)];
    else         cur = &lft[n_pad-1 + n_pad*k];
  }
  else
  {
    if (k>=kmax) cur = &top[k-kmax + n_pad*(i+istart)];
    else         cur = &rings[i + iblk*k];
  }
  jj = j_bound(j - (*cur).rot_j);

  ut = (*cur).C[jj].u;
  vt = (*cur).C[jj].v;
  wt = (*cur).C[jj].w;
  dz = (*cur).zc;
  k--;

  ///////////////////////////////////////

  k--;
  if (i==iblk) 
  {
    if (k==-1) cur = &rgh[0 + n_pad*(k+1)];
    else       cur = &rgh[0 + n_pad*k];
  }
  else if (i==-1)
  {
    if (k==-1) cur = &lft[n_pad-1 + n_pad*(k+1)];
    else       cur = &lft[n_pad-1 + n_pad*k];
  }
  else
  {
    if (k<=-1) cur = &udr[k+n_pad + n_pad*(i+istart)];
    else       cur = &rings[i + iblk*k];
  }
  jj = j_bound(j - (*cur).rot_j);

  ud = (*cur).C[jj].u;
  vd = (*cur).C[jj].v;
  wd = (*cur).C[jj].w;
  dz-= (*cur).zc;

  #endif

  return;
}

__global__ void viscous_tensor(double* tau, int istart, int iblk, ring_geometry *ygrid, ring *rings, ring *lft, ring *rgh, ring *udr, ring *top)
{
	int n = threadIdx.x + blockDim.x*blockIdx.x;

	int i = n/jmax;
	int k = blockIdx.y;

	if (n>=(iblk+2)*jmax) return;
	#if ndim==3
	if (i==0 && k==0) return;
	if (i==iblk+1 && k==kmax+1) return;
	if (i==0 && k==kmax+1) return;
	if (i==iblk+1 && k==0) return;
	#endif

	double r;
	double dr, dphi, dz;
	double rho, u, v, w;
	double uf, ub, ur, ul, ut, ud;
	double vf, vb, vr, vl, vt, vd;
	double wf, wb, wr, wl, wt, wd;

	load_cross(istart, iblk, ygrid, rings, lft, rgh, udr, top,
	           r, dr, dphi, dz, rho, u, v, w, 
	           ur, vr, wr, ul, vl, wl, uf, vf, wf,
	           ub, vb, wb, ut, vt, wt, ud, vd, wd);

	double nr = get_nu(r)*rho;

	double dudx = (uf-ub)/dr;
	double dvdx = (vf-vb)/dr;
	double dwdx = (wf-wb)/dr;

	double dudy = (ur-ul)/dphi/r;
	double dvdy = (vr-vl)/dphi/r;
	double dwdy = (wr-wl)/dphi/r;
  
	#if ndim>2
	double dudz = (ut-ud)/dz/r;
	double dvdz = (vt-vd)/dz/r;
	double dwdz = (wt-wd)/dz/r;
	#endif

	int ind = (n + k*(iblk+2)*jmax)*6;

	tau[ind+0] = 2.0*nr*(dudx);
	tau[ind+1] =     nr*(dudy + dvdx - v/r);
	tau[ind+3] = 2.0*nr*(dvdy + u/r);
	#if ndim>2
	tau[ind+2] =     nr*(dudz + dwdx);
	tau[ind+4] =     nr*(dvdz + dwdy);
	tau[ind+5] = 2.0*nr*(dwdz);
	#endif

	return;
}

__global__ void viscous_force(double* tau, int iblk, ring_geometry *ygrid, ring *rings, double dt)
{
	int n = threadIdx.x + blockDim.x*blockIdx.x;

	int i = n/jmax;
	int j = n - i*jmax;
	int k = blockIdx.y;

	#if ndim>2
	int ind = (j + (i+1)*jmax + (k+1)*(iblk+2)*jmax)*6;
	#else
	int ind = (j + (i+1)*jmax)*6;
	#endif

	double* t0 = &tau[ind];
	double* tf = &tau[ind+jmax*6];
	double* tb = &tau[ind-jmax*6];
	double* tr = &tau[ind+6];
	double* tl = &tau[ind-6];
	#if ndim>2
	double* tt = &tau[ind+(iblk+2)*jmax*6];
	double* td = &tau[ind-(iblk+2)*jmax*6];
	#endif

	double rad, pol, rad_cyl;
	double dr, dphi, dz;
	double fx, fy, fz;

	int jj;
	int ik = i + iblk*k;
	Cell C_tmp;

	if (i<iblk && j<jmax)
	{
		jj = j_bound(j - rings[ik].rot_j);
		assign_pos(rad, pol, rad_cyl, rings[ik].xc, rings[ik].zc);
		dr   = 2.0*rings[ik].dx;
		dphi = 2.0*(*ygrid).dy[j]*rad_cyl;
		dz   = 2.0*rings[ik].dz*rad;

		C_tmp = rings[ik].C[jj];

		fx = t0[0]/rad + (tf[0]-tb[0])/dr + (tr[1]-tl[1])/dphi - t0[3]/rad;
		fy = 2.0*t0[1]/rad + (tf[1]-tb[1])/dr + (tr[3]-tl[3])/dphi;
		fz = 0.0;

		#if ndim>2
		fx += (tt[2]-td[2])/dz;
		fy += (tt[4]-td[4])/dz;
		fz += t0[2]/rad + (tf[2]-tb[2])/dr + (tr[4]-tl[4])/dphi + (tt[5]-td[5])/dz;
		#endif

		C_tmp.u += fx*dt/C_tmp.r;
		C_tmp.v += fy*dt/C_tmp.r;
		C_tmp.w += fz*dt/C_tmp.r;

		rings[ik].C[jj] = C_tmp;
	}
	return;
}
