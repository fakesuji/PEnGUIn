__device__ void set_state(int i, int geom, double* xa, double* dx, double* dv, double rad, 
                          double* r, double* p, double* u, double* v, double* w, double dt, double us, State &S,
                          bool print=false)
{
	double rl_par[4], pl_par[4], ul_par[4];
	double rr_par[4], pr_par[4], ur_par[4];

	double rl, pl, ul, cl;
	double rr, pr, ur, cr;
	double xl, xr, tmp, dis;
	double r0, p0, u0;
	double p_, u_;
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

	double rmin, pmin;
	rmin = r[i]*1.0e-10;
	pmin = p[i]*1.0e-10;

	//rr = get_CON_aveL(geom, xl, tmp, xr, rr_par);
	//pr = get_CON_aveL(geom, xl, tmp, xr, pr_par);
	ur = get_PRM_aveL(geom, xl, tmp, xr, ur_par);
	//rr = fmax(rr,rmin);
	//pr = fmax(pr,pmin);

	dis = fmin(ur, 0.0)*dt;
	//cr = sqrt(gam*pr/rr); 

	/////////////////////////////////////////////////////////

	tmp = xl - dis;
	if (tmp>xr) tmp = xr;

	r0 = get_CON_aveL(geom, xl, tmp, xr, rr_par);
	p0 = get_CON_aveL(geom, xl, tmp, xr, pr_par);
	u0 = get_PRM_aveL(geom, xl, tmp, xr, ur_par);
	r0 = fmax(r0, rmin);
	p0 = fmax(p0, pmin);
	cr = sqrt(gam*p0/r0); 

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

		u_ = 0.5*(up-um)/cr;
		p_ = 0.5*(pp-pm)/r0/cr;
	}
	else
	{
		u_ = (u0-um)/cr;
		p_ = (p0-pm)/r0/cr;
	}
	
	S.rr = r0*exp_lim(u_);
	S.pr = p0*exp_lim(gam*u_);
	S.ur = u0 + p_ + us;

	if (print) printf(" rr=%e, pr=%e, ur=%e\n r0=%e, p0%e, u0=%e\n u_=%e, p_=%e",rr,pr,ur,r0,p0,u0,u_,p_);

	/////////////////////////////////////////////////////////

	get_PRM_parameters(i, geom, xa, dx, dv, v, rr_par);
	get_PRM_parameters(i, geom, xa, dx, dv, w, pr_par);
	tmp = xl - dis;
	if (tmp>xr) tmp = xr;
	S.vr = get_PRM_aveL(geom, xl, tmp, xr, rr_par);
	S.wr = get_PRM_aveL(geom, xl, tmp, xr, pr_par);

	/////////////////////////////////////////////////////////

	xl = xa[i-1];
	xr = xa[i];

	tmp = fmax(u[i-1],0.0);
	tmp = xr - (sqrt(gam*p[i-1]/r[i-1]) + tmp)*dt;
	if (tmp<xl) tmp = xl;

	rmin = r[i-1]*1.0e-10;
	pmin = p[i-1]*1.0e-10;

	//rl = get_CON_aveR(geom, xl, tmp, xr, rl_par);
	//pl = get_CON_aveR(geom, xl, tmp, xr, pl_par);
	ul = get_PRM_aveR(geom, xl, tmp, xr, ul_par);
	//rl = fmax(rl, rmin);
	//pl = fmax(pl, pmin);

	dis = fmax(ul,0.0)*dt;
	//cl = sqrt(gam*pl/rl);

	/////////////////////////////////////////////////////////

	tmp = xr - dis;
	if (tmp<xl) tmp = xl;

	r0 = get_CON_aveR(geom, xl, tmp, xr, rl_par);
	p0 = get_CON_aveR(geom, xl, tmp, xr, pl_par);
	u0 = get_PRM_aveR(geom, xl, tmp, xr, ul_par);
	r0 = fmax(r0, rmin);
	p0 = fmax(p0, pmin);
	cl = sqrt(gam*p0/r0);

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

		u_ = 0.5*(up-um)/cl;
		p_ = 0.5*(pp-pm)/r0/cl;
	}
	else
	{
		u_ = (up-u0)/cl;
		p_ = (pp-p0)/r0/cl;
	}

	S.rl = r0*exp_lim(u_);
	S.pl = p0*exp_lim(gam*u_);
	S.ul = u0 + p_ + us;

	if (print) printf(" rl=%e, pl=%e, ul=%e\n r0=%e, p0%e, u0=%e\n u_=%e, p_=%e",rl,pl,ul,r0,p0,u0,u_,p_);

	/////////////////////////////////////////////////////////

	get_PRM_parameters(i-1, geom, xa, dx, dv, v, rl_par);
	get_PRM_parameters(i-1, geom, xa, dx, dv, w, pl_par);
	tmp = xr - dis;
	if (tmp<xl) tmp = xl;
	S.vl = get_PRM_aveR(geom, xl, tmp, xr, rl_par);
	S.wl = get_PRM_aveR(geom, xl, tmp, xr, pl_par);

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
