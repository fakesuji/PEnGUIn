__device__ void set_state(int i, int geom, double* xa, double* dx, double* dv, double rad, 
                          double* r, double* p, double* u, double* v, double* w, double dt, double us, State &S,
                          bool print=false)
{
	double r_par[4], p_par[4], u_par[4];

	double ul, ur;
	double xl, xr, tmp, dis, q;
	double ql, qr;
	double r0, p0, u0;
	double p_, u_;

	/////////////////////////////////////////////////////////

	get_CON_parameters(i, geom, xa, dx, dv, r, r_par);
	get_CON_parameters(i, geom, xa, dx, dv, p, p_par);
	get_PRM_parameters(i, geom, xa, dx, dv, u, u_par);

	xl = xa[i];
	xr = xa[i+1];

	tmp = fmin(u_par[1]-u_par[0],0.0);
	tmp = xl - tmp*dt + sqrt(gam*p[i]/r[i])*dt;
	dimensionless_x(xl,tmp,xr,q,ql,qr);

	double rmin, pmin;
	rmin = r[i]*1.0e-10;
	pmin = p[i]*1.0e-10;

	ur = get_PRM_aveL(geom, q, u_par, ql, qr);
	dis = fmin(ur, 0.0)*dt;

	/////////////////////////////////////////////////////////

	tmp = xl - dis;
	dimensionless_x(xl,tmp,xr,q,ql,qr);

	r0 = get_CON_aveL(geom, q, r_par, ql, qr);
	p0 = get_CON_aveL(geom, q, p_par, ql, qr);
	u0 = get_PRM_aveL(geom, q, u_par, ql, qr);
	r0 = fmax(r0, rmin);
	p0 = fmax(p0, pmin);

	u_ = 0.5*dt*(u0-u[i])/(0.5*(xr-tmp));
	p_ = 0.5*dt*(p0-p[i])/(0.5*(xr-tmp))/r0;
	
	S.rr = r0*exp_lim(u_);
	S.pr = p0*exp_lim(gam*u_);
	S.ur = u0 + p_ + us;

	if (print) printf(" ur=%e\n r0=%e, p0=%e, u0=%e\n u_=%e, p_=%e\n",ur,r0,p0,u0,u_,p_);

	/////////////////////////////////////////////////////////

	get_PRM_parameters(i, geom, xa, dx, dv, v, r_par);
	get_PRM_parameters(i, geom, xa, dx, dv, w, p_par);

	S.vr = get_PRM_aveL(geom, q, r_par, ql, qr);
	S.wr = get_PRM_aveL(geom, q, p_par, ql, qr);

	/////////////////////////////////////////////////////////

	get_CON_parameters(i-1, geom, xa, dx, dv, r, r_par);
	get_CON_parameters(i-1, geom, xa, dx, dv, p, p_par);
	get_PRM_parameters(i-1, geom, xa, dx, dv, u, u_par);

	xl = xa[i-1];
	xr = xa[i];

	tmp = fmax(u_par[1]+u_par[2],0.0);
	tmp = xr - tmp*dt - sqrt(gam*p[i-1]/r[i-1])*dt;
	dimensionless_x(xl,tmp,xr,q,ql,qr);

	rmin = r[i-1]*1.0e-10;
	pmin = p[i-1]*1.0e-10;

	ul = get_PRM_aveR(geom, q, u_par, ql, qr);
	dis = fmax(ul,0.0)*dt;

	/////////////////////////////////////////////////////////

	tmp = xr - dis;
	dimensionless_x(xl,tmp,xr,q,ql,qr);

	r0 = get_CON_aveR(geom, q, r_par, ql, qr);
	p0 = get_CON_aveR(geom, q, p_par, ql, qr);
	u0 = get_PRM_aveR(geom, q, u_par, ql, qr);
	r0 = fmax(r0, rmin);
	p0 = fmax(p0, pmin);

	u_ = 0.5*dt*(u[i-1]-u0)/(0.5*(tmp-xl));
	p_ = 0.5*dt*(p[i-1]-p0)/(0.5*(tmp-xl))/r0;

	S.rl = r0*exp_lim(u_);
	S.pl = p0*exp_lim(gam*u_);
	S.ul = u0 + p_ + us;

	if (print) printf(" ul=%e\n r0=%e, p0=%e, u0=%e\n u_=%e, p_=%e\n",ul,r0,p0,u0,u_,p_);

	/////////////////////////////////////////////////////////

	get_PRM_parameters(i-1, geom, xa, dx, dv, v, r_par);
	get_PRM_parameters(i-1, geom, xa, dx, dv, w, p_par);

	S.vl = get_PRM_aveR(geom, q, r_par, ql, qr);
	S.wl = get_PRM_aveR(geom, q, p_par, ql, qr);

	return;
}
