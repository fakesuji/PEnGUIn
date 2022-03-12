__device__ void set_L_state(int i, int geom, double* xa, double* dx, double* dv, double rad, 
                            double* r, double* p, double* u, double* v, double* w, double dt, double force, State &S)
{
	double r_par[4], p_par[4], u_par[4];
	get_CON_parameters(i, geom, xa, dx, dv, r, r_par);
	get_CON_parameters(i, geom, xa, dx, dv, p, p_par);
	get_PRM_parameters(i, geom, xa, dx, dv, u, u_par);

	double xl = xa[i];
	double xr = xa[i+1];
	double cl = sqrt(gam*p[i]/r[i]);
	double tmp = xr - fmax(cl, u[i]+cl)*dt;

	double rl = get_CON_aveR(geom, xl, tmp, xr, r_par);
	double pl = get_CON_aveR(geom, xl, tmp, xr, p_par);
	double ul = get_PRM_aveR(geom, xl, tmp, xr, u_par);
	double C = sqrt(gam*pl*rl);
	cl = C/rl;

	double r_, p_, u_;

	double Bp = 0.0;
	double Bm = 0.0;
	double B0 = 0.0;

	if (ul>-cl)
	{
		tmp = xr - (ul+cl)*dt;
		p_ = get_CON_aveR(geom, xl, tmp, xr, p_par);
		u_ = get_PRM_aveR(geom, xl, tmp, xr, u_par);

		Bp = -( (ul-u_) + (pl-p_)/C - dt*force)/2.0/C;
	}

	if (ul>0.0)
	{
		tmp = xr - ul*dt;
		r_ = get_CON_aveR(geom, xl, tmp, xr, r_par);
		p_ = get_CON_aveR(geom, xl, tmp, xr, p_par);

		B0 = (pl-p_)/C/C + 1.0/rl - 1.0/r_;
	}

	if (ul>cl)
	{
		tmp = xr - (ul-cl)*dt;
		p_ = get_CON_aveR(geom, xl, tmp, xr, p_par);
		u_ = get_PRM_aveR(geom, xl, tmp, xr, u_par);

		Bm =  ( (ul-u_) - (pl-p_)/C - dt*force)/2.0/C;
	}

	S.pl = pl + (Bp + Bm)*C*C;
	S.ul = ul + (Bp - Bm)*C;
	S.rl = 1.0/( 1.0/rl - (B0 + Bp + Bm) );

	return;
}

__device__ void set_R_state(int i, int geom, double* xa, double* dx, double* dv, double rad, 
                            double* r, double* p, double* u, double* v, double* w, double dt, double force, State &S)
{
	double r_par[4], p_par[4], u_par[4];
	get_CON_parameters(i, geom, xa, dx, dv, r, r_par);
	get_CON_parameters(i, geom, xa, dx, dv, p, p_par);
	get_PRM_parameters(i, geom, xa, dx, dv, u, u_par);

	double xl = xa[i];
	double xr = xa[i+1];
	double cr = sqrt(gam*p[i]/r[i]);
	double tmp = xl + fmax(cr, -u[i]+cr)*dt;

	double rr = get_CON_aveL(geom, xl, tmp, xr, r_par);
	double pr = get_CON_aveL(geom, xl, tmp, xr, p_par);
	double ur = get_PRM_aveL(geom, xl, tmp, xr, u_par);
	double C = sqrt(gam*pr*rr);
	cr = C/rr;

	double r_, p_, u_;

	double Bp = 0.0;
	double Bm = 0.0;
	double B0 = 0.0;

	if (ur<-cr)
	{
		tmp = xl - (ur+cr)*dt;
		p_ = get_CON_aveL(geom, xl, tmp, xr, p_par);
		u_ = get_PRM_aveL(geom, xl, tmp, xr, u_par);

		Bp = -( (ur-u_) + (pr-p_)/C - dt*force)/2.0/C;
	}

	if (ur<0.0)
	{
		tmp = xl - ur*dt;
		r_ = get_CON_aveL(geom, xl, tmp, xr, r_par);
		p_ = get_CON_aveL(geom, xl, tmp, xr, p_par);

		B0 = (pr-p_)/C/C + 1.0/rr - 1.0/r_;
	}

	if (ur<cr)
	{
		tmp = xl - (ur-cr)*dt;
		p_ = get_CON_aveL(geom, xl, tmp, xr, p_par);
		u_ = get_PRM_aveL(geom, xl, tmp, xr, u_par);

		Bm =  ( (ur-u_) - (pr-p_)/C - dt*force)/2.0/C;
	}

	S.pr = pr + (Bp + Bm)*C*C;
	S.ur = ur + (Bp - Bm)*C;
	S.rr = 1.0/( 1.0/rr - (B0 + Bp + Bm) );

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
