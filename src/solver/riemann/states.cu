
__device__ void set_state(int i, int geom, double* xa, double* dx, double* dv, double rad, 
                          double* r, double* p, double* u, double* v, double* w, double dt, double us, State &S)
{
	double rl_par[4], pl_par[4], ul_par[4];
	double rr_par[4], pr_par[4], ur_par[4];

	double rl, pl, ul, cl;
	double rr, pr, ur, cr;
	double xl, xr, tmp, dis;
	double r_, p_, u_;
	double r0, p0, u0, c0;
	double pp, up;
	double pm, um;

	get_CON_parameters(i-1, geom, xa, dx, dv, r, rl_par);
	get_CON_parameters(i-1, geom, xa, dx, dv, p, pl_par);
	get_PRM_parameters(i-1, geom, xa, dx, dv, u, ul_par);

	get_CON_parameters(i, geom, xa, dx, dv, r, rr_par);
	get_CON_parameters(i, geom, xa, dx, dv, p, pr_par);
	get_PRM_parameters(i, geom, xa, dx, dv, u, ur_par);

	/////////////////////////////////////////////////////////

	xl = xa[i];
	xr = xa[i+1];

	tmp = fmax(-u[i],0.0);
	tmp = xl + (sqrt(gam*p[i]/r[i]) + tmp)*dt;
	if (tmp>xr) tmp = xr;

	rr = get_CON_aveL(geom, xl, tmp, xr, rr_par);
	pr = get_CON_aveL(geom, xl, tmp, xr, pr_par);
	ur = get_PRM_aveL(geom, xl, tmp, xr, ur_par);

	dis = fmin(ur, 0.0)*dt;
	cr = sqrt(gam*pr/rr);

	/////////////////////////////////////////////////////////

	tmp = xl - dis;
	if (tmp>xr) tmp = xr;

	r0 = get_CON_aveL(geom, xl, tmp, xr, rr_par);
	p0 = get_CON_aveL(geom, xl, tmp, xr, pr_par);
	u0 = get_PRM_aveL(geom, xl, tmp, xr, ur_par);
	c0 = sqrt(gam*p0/r0);

	tmp = xl - dis + cr*dt;
	if (tmp>xr) tmp = xr;

	pm = get_CON_aveL(geom, xl, tmp, xr, pr_par);
	um = get_PRM_aveL(geom, xl, tmp, xr, ur_par);

	tmp = xl - dis - cr*dt;
	if (tmp>xr) tmp = xr;

	if (tmp>xl)
	{
		pp = get_CON_aveL(geom, xl, tmp, xr, pr_par);
		up = get_PRM_aveL(geom, xl, tmp, xr, ur_par);

		S.rr = r0 + 0.5*(up-um)*r0/c0;// - 0.5*(pp-pm)/c0/c0;
		S.pr = p0 + 0.5*(up-um)*r0*c0;
		S.ur = u0 + 0.5*(pp-pm)/r0/c0 + 0.5*us;
	}
	else
	{
		S.rr = r0 + (u0-um)*r0/c0;// - (p0-pm)/c0/c0;
		S.pr = p0 + (u0-um)*r0*c0;
		S.ur = u0 + (p0-pm)/r0/c0 + 0.5*us;
	}

	/////////////////////////////////////////////////////////

	xl = xa[i-1];
	xr = xa[i];

	tmp = fmax(u[i-1],0.0);
	tmp = xr - (sqrt(gam*p[i-1]/r[i-1]) + tmp)*dt;
	if (tmp<xl) tmp = xl;

	rl = get_CON_aveR(geom, xl, tmp, xr, rl_par);
	pl = get_CON_aveR(geom, xl, tmp, xr, pl_par);
	ul = get_PRM_aveR(geom, xl, tmp, xr, ul_par);
	dis = fmax(ul,0.0)*dt;
	cl = sqrt(gam*pl/rl);

	/////////////////////////////////////////////////////////

	tmp = xr - dis;
	if (tmp<xl) tmp = xl;

	r0 = get_CON_aveR(geom, xl, tmp, xr, rl_par);
	p0 = get_CON_aveR(geom, xl, tmp, xr, pl_par);
	u0 = get_PRM_aveR(geom, xl, tmp, xr, ul_par);
	c0 = sqrt(gam*p0/r0);

	tmp = xr - dis - cl*dt;
	if (tmp<xl) tmp = xl;

	pp = get_CON_aveR(geom, xl, tmp, xr, pl_par);
	up = get_PRM_aveR(geom, xl, tmp, xr, ul_par);

	tmp = xr - dis + cl*dt;
	if (tmp<xl) tmp = xl;

	if (tmp<xr)
	{
		pm = get_CON_aveR(geom, xl, tmp, xr, pl_par);
		um = get_PRM_aveR(geom, xl, tmp, xr, ul_par);

		S.rl = r0 + 0.5*(up-um)*r0/c0;// + 0.5*(pp-pm)/c0/c0;
		S.pl = p0 + 0.5*(up-um)*r0*c0;
		S.ul = u0 + 0.5*(pp-pm)/r0/c0 + 0.5*us;
	}
	else
	{
		S.rl = r0 + (up-u0)*r0/c0;// + (pp-p0)/c0/c0;
		S.pl = p0 + (up-u0)*r0*c0;
		S.ul = u0 + (pp-p0)/r0/c0 + 0.5*us;
	}

	

	S.pl = fmax(smallp*p[i-1],S.pl);
	S.pr = fmax(smallp*p[i],S.pr);

	S.rl = fmax(smallr*r[i-1],S.rl);
	S.rr = fmax(smallr*r[i],S.rr);

	return;
}
/*
__device__ void set_state(int i, int geom, double* xa, double* dx, double* dv, double rad, 
                          double* r, double* p, double* u, double* v, double* w, double dt, double us, State &S)
{
	double rl_par[4], pl_par[4], ul_par[4];
	double rr_par[4], pr_par[4], ur_par[4];

	double rl, pl, ul, cl;
	double rr, pr, ur, cr;
	double xl, xr, tmp;
	double r_, p_, u_;
	double Bp, Bm, B0;

	get_CON_parameters(i-1, geom, xa, dx, dv, r, rl_par);
	get_CON_parameters(i-1, geom, xa, dx, dv, p, pl_par);
	get_PRM_parameters(i-1, geom, xa, dx, dv, u, ul_par);

	get_CON_parameters(i, geom, xa, dx, dv, r, rr_par);
	get_CON_parameters(i, geom, xa, dx, dv, p, pr_par);
	get_PRM_parameters(i, geom, xa, dx, dv, u, ur_par);

	/////////////////////////////////////////////////////////

	xl = xa[i-1];
	xr = xa[i];

	tmp = fmax(u[i-1],0.0);
	tmp = xr - (sqrt(gam*p[i-1]/r[i-1]) + tmp)*dt;

	if (tmp<xl) 
	{
		tmp = xl;
	}

	rl = get_CON_aveR(geom, xl, tmp, xr, rl_par);
	pl = get_CON_aveR(geom, xl, tmp, xr, pl_par);
	ul = get_PRM_aveR(geom, xl, tmp, xr, ul_par);
	cl = sqrt(gam*pl/rl);

	Bp = 0.0;
	Bm = 0.0;
	B0 = 0.0;

	if (ul>-cl)
	{
		tmp = xr - (ul+cl)*dt;

		p_ = get_CON_aveR(geom, xl, tmp, xr, pl_par);
		u_ = get_PRM_aveR(geom, xl, tmp, xr, ul_par);

		Bp = -( (ul-u_)/cl + (pl-p_)/pl - us/cl)/2.0;
	}

	if (ul>0.0)
	{
		tmp = xr - ul*dt;

		r_ = get_CON_aveR(geom, xl, tmp, xr, rl_par);
		p_ = get_CON_aveR(geom, xl, tmp, xr, pl_par);

		B0 = 2.0 - rl/r_ - p_/pl;
	}

	if (ul>cl)
	{
		tmp = xr - (ul-cl)*dt;

		p_ = get_CON_aveR(geom, xl, tmp, xr, pl_par);
		u_ = get_PRM_aveR(geom, xl, tmp, xr, ul_par);

		Bm =  ( (ul-u_)/cl - (pl-p_)/pl - us/cl)/2.0;
	}

	S.rl = rl / (1.0 - B0 - Bp - Bm);
	S.pl = pl * (1.0 + Bp + Bm);
	S.ul = ul + Bp*cl - Bm*cl;

	/////////////////////////////////////////////////////////

	xl = xa[i];
	xr = xa[i+1];

	tmp = fmax(-u[i],0.0);
	tmp = xl + (sqrt(gam*p[i]/r[i]) + tmp)*dt;

	if (tmp>xr) 
	{
		tmp = xr;
	}

	rr = get_CON_aveL(geom, xl, tmp, xr, rr_par);
	pr = get_CON_aveL(geom, xl, tmp, xr, pr_par);
	ur = get_PRM_aveL(geom, xl, tmp, xr, ur_par);
	cr = sqrt(gam*pr/rr);

	Bp = 0.0;
	Bm = 0.0;
	B0 = 0.0;

	if (ur<-cr)
	{
		tmp = xl - (ur+cr)*dt;

		p_ = get_CON_aveL(geom, xl, tmp, xr, pr_par);
		u_ = get_PRM_aveL(geom, xl, tmp, xr, ur_par);

		Bp = -( (ur-u_)/cr + (pr-p_)/pr - us/cr)/2.0;
	}

	if (ur<0.0)
	{
		tmp = xl - ur*dt;

		r_ = get_CON_aveL(geom, xl, tmp, xr, rr_par);
		p_ = get_CON_aveL(geom, xl, tmp, xr, pr_par);

		B0 = 2.0 - rr/r_ - p_/pr;
	}

	if (ur<cr)
	{
		tmp = xl - (ur-cr)*dt;

		p_ = get_CON_aveL(geom, xl, tmp, xr, pr_par);
		u_ = get_PRM_aveL(geom, xl, tmp, xr, ur_par);

		Bm =  ( (ur-u_)/cr - (pr-p_)/pr - us/cr)/2.0;
	}

	S.rr = rr / (1.0 - B0 - Bp - Bm);
	S.pr = pr * (1.0 + Bp + Bm);
	S.ur = ur + Bp*cr - Bm*cr;

	S.pl = fmax(smallp*p[i-1],S.pl);
	S.rl = fmax(smallr*r[i-1],S.rl);
	S.pr = fmax(smallp*p[i],S.pr);
	S.rr = fmax(smallr*r[i],S.rr);

	return;
}
*/
__device__ void set_R_state_passive(int i, int geom, double* xa, double* dx, double* dv, double rad, 
                                    double sr, double sm, double* u, double* v, double* w, double dt, State &S)
{
	double v_par[4], w_par[4];
	get_PRM_parameters(i, geom, xa, dx, dv, v, v_par);
	get_PRM_parameters(i, geom, xa, dx, dv, w, w_par);

	double tmp;
	double xl = xa[i];
	double xr = xa[i+1];

	if (signbit(sr)!=0)
	{
		tmp = xl - u[i]*dt;
		S.vr = get_PRM_aveL(geom, xl, tmp, xr, v_par);
		S.wr = get_PRM_aveL(geom, xl, tmp, xr, w_par);
	}
	else if (signbit(sm)!=0)
	{
		tmp = xl - sm*dt;
		S.vr = get_PRM_aveL(geom, xl, tmp, xr, v_par);
		S.wr = get_PRM_aveL(geom, xl, tmp, xr, w_par);
	}
	else
	{
		S.vr = 0.0;
		S.wr = 0.0;
	}

	return;
}

__device__ void set_L_state_passive(int i, int geom, double* xa, double* dx, double* dv, double rad, 
                                    double sl, double sm, double* u, double* v, double* w, double dt, State &S)
{
	double v_par[4], w_par[4];
	get_PRM_parameters(i, geom, xa, dx, dv, v, v_par);
	get_PRM_parameters(i, geom, xa, dx, dv, w, w_par);

	double tmp;
	double xl = xa[i];
	double xr = xa[i+1];

	if (signbit(sl)==0)
	{
		tmp = xr - u[i]*dt;
		S.vl = get_PRM_aveR(geom, xl, tmp, xr, v_par);
		S.wl = get_PRM_aveR(geom, xl, tmp, xr, w_par);
	}
	else if (signbit(sm)==0)
	{
		tmp = xr - sm*dt;
		S.vl = get_PRM_aveR(geom, xl, tmp, xr, v_par);
		S.wl = get_PRM_aveR(geom, xl, tmp, xr, w_par);
	}
	else
	{
		S.vl = 0.0;
		S.wl = 0.0;
	}

	return;
}
