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

	#if recon_flag==0
	if ( rl_par[2]!=rl_par[1] && rr_par[0]!=rr_par[1] )
	{
		tmp = 0.5*(rl_par[2]+rr_par[0]);
		rl_par[2] = tmp;
		rr_par[0] = tmp;
	}
	if ( pl_par[2]!=pl_par[1] && pr_par[0]!=pr_par[1] )
	{
		tmp = 0.5*(pl_par[2]+pr_par[0]);
		pl_par[2] = tmp;
		pr_par[0] = tmp;
	}
	#endif

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

		//B0 = (1.0/rl/gam-p_/C/C) + 1.0/rl - 1.0/r_;
		B0 = 2.0 - rl/r_ - p_/pl;
	}

	if (ul>cl)
	{
		tmp = xr - (ul-cl)*dt;

		p_ = get_CON_aveR(geom, xl, tmp, xr, pl_par);
		u_ = get_PRM_aveR(geom, xl, tmp, xr, ul_par);

		Bm =  ( (ul-u_)/cl - (pl-p_)/pl - us/cl)/2.0;
	}

	S.rl = rl/( 1.0 - B0 - Bp - Bm );
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

		//B0 = (1.0/rr/gam-p_/C/C) + 1.0/rr - 1.0/r_;
		B0 = 2.0 - rr/r_ - p_/pr;
	}

	if (ur<cr)
	{
		tmp = xl - (ur-cr)*dt;

		p_ = get_CON_aveL(geom, xl, tmp, xr, pr_par);
		u_ = get_PRM_aveL(geom, xl, tmp, xr, ur_par);

		Bm =  ( (ur-u_)/cr - (pr-p_)/pr - us/cr)/2.0;
	}

	S.rr = rr/( 1.0 - B0 - Bp - Bm );
	S.pr = pr * (1.0 + Bp + Bm);
	S.ur = ur + Bp*cr - Bm*cr;

	S.pl = fmax(smallp*p[i-1],S.pl);
	S.rl = fmax(smallr*r[i-1],S.rl);
	S.pr = fmax(smallp*p[i],S.pr);
	S.rr = fmax(smallr*r[i],S.rr);

	return;
}

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
