__device__ void set_state(int i, int geom, double* xa, double* dx, double* dv, double rad, 
                          double* r, double* p, double* u, double* v, double* w, double dt, double us, State &S,
                          bool print=false)
{
	double r_par[4], p_par[4], u_par[4];

	double ul, ur;
	double xl, xr, tmp, dis, q;
	double ql, qr;
	double r0, p0, u0;
	double p_, u_, cs;
	double hdt = 0.5*dt;

	/////////////////////////////////////////////////////////

	get_CON_parameters(i, geom, xa, dx, dv, r, r_par);
	get_CON_parameters(i, geom, xa, dx, dv, p, p_par);
	get_PRM_parameters(i, geom, xa, dx, dv, u, u_par);

	xl = xa[i];
	xr = xa[i+1];

	cs = sqrt(gam*p[i]/r[i]);
	tmp = fmin(u[i],0.0);
	tmp = xl - tmp*dt + cs*dt;
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
	cs = sqrt(gam*p0/r0);

	tmp = xl - dis + cs*dt;
	dimensionless_x(xl,tmp,xr,q,ql,qr);

	u_ = get_PRM_aveL(geom, q, u_par, ql, qr);
	p_ = get_CON_aveL(geom, q, p_par, ql, qr);

	u_ = 2.0*(u_-u0)/(cs*dt);
	p_ = 2.0*(p_-p0)/(cs*dt);

	u_ *= -hdt;
	p_ *= -hdt;
	
	S.rr = r0*exp_lim(u_);
	S.pr = p0*exp_lim(gam*u_);
	S.ur = u0 + p_/r0 + us;

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

	cs = sqrt(gam*p[i-1]/r[i-1]);
	tmp = fmax(u[i-1],0.0);
	tmp = xr - tmp*dt - cs*dt;
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
	cs = sqrt(gam*p0/r0);

	tmp = xr - dis - cs*dt;
	dimensionless_x(xl,tmp,xr,q,ql,qr);

	u_ = get_PRM_aveR(geom, q, u_par, ql, qr);
	p_ = get_CON_aveR(geom, q, p_par, ql, qr);

	u_ = 2.0*(u0-u_)/(cs*dt);
	p_ = 2.0*(p0-p_)/(cs*dt);

	u_ *= -hdt;
	p_ *= -hdt;

	S.rl = r0*exp_lim(u_);
	S.pl = p0*exp_lim(gam*u_);
	S.ul = u0 + p_/r0 + us;

	if (print) printf(" ul=%e\n r0=%e, p0=%e, u0=%e\n u_=%e, p_=%e\n",ul,r0,p0,u0,u_,p_);

	/////////////////////////////////////////////////////////

	get_PRM_parameters(i-1, geom, xa, dx, dv, v, r_par);
	get_PRM_parameters(i-1, geom, xa, dx, dv, w, p_par);

	S.vl = get_PRM_aveR(geom, q, r_par, ql, qr);
	S.wl = get_PRM_aveR(geom, q, p_par, ql, qr);

	return;
}
