__device__ void set_L_state(int i, int geom, double* xa, double* dx, double* dv, double rad, 
                            double* r, double* p, double* u, double* v, double* w, double dt, State &S)
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

	double B_ = 0.0;

	S.rl = rl;
	S.pl = pl;
	S.ul = ul;
//if (xa[i]<2.0 && xa[i]>1.993) printf("%f, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e\n", rad, dx[i-1], dx[i], dx[i+1], u[i-1], u_par[0], u_par[1], u_par[2], u[i+1]);
	if (ul>cl)
	{
		tmp = xr - (ul-cl)*dt;
		p_ = pl - get_CON_aveR(geom, xl, tmp, xr, p_par);
		u_ = ul - get_PRM_aveR(geom, xl, tmp, xr, u_par);

		B_ =  (u_*rl*cl - p_)/2.0;

		S.rl -= -B_/(cl*cl);
		S.pl -= -B_;
		S.ul -=  B_/(rl*cl);
	}

	if (ul>0.0)
	{
		tmp = xr - ul*dt;
		r_ = rl - get_CON_aveR(geom, xl, tmp, xr, r_par);
		p_ = pl - get_CON_aveR(geom, xl, tmp, xr, p_par);

		S.rl -= r_ - p_/(cl*cl);
	}

	if (ul>-cl)
	{
		tmp = xr - (ul+cl)*dt;
		p_ = pl - get_CON_aveR(geom, xl, tmp, xr, p_par);
		u_ = ul - get_PRM_aveR(geom, xl, tmp, xr, u_par);

		B_ = (u_*rl*cl + p_)/2.0;

		S.rl -= B_/(cl*cl);
		S.pl -= B_;
		S.ul -= B_/(rl*cl);
	}

	return;
}

__device__ void set_R_state(int i, int geom, double* xa, double* dx, double* dv, double rad, 
                            double* r, double* p, double* u, double* v, double* w, double dt, State &S)
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

	double B_ = 0.0;
	
	S.rr = rr;
	S.pr = pr;
	S.ur = ur;

	if (ur<-cr)
	{
		tmp = xl - (ur+cr)*dt;
		p_ = pr - get_CON_aveL(geom, xl, tmp, xr, p_par);
		u_ = ur - get_PRM_aveL(geom, xl, tmp, xr, u_par);

		B_ = (u_*rr*cr + p_)/2.0;
	
		S.rr -= B_/(cr*cr);
		S.pr -= B_;
		S.ur -= B_/(rr*cr);
	}

	if (ur<0.0)
	{
		tmp = xl - ur*dt;
		r_ = rr - get_CON_aveL(geom, xl, tmp, xr, r_par);
		p_ = pr - get_CON_aveL(geom, xl, tmp, xr, p_par);

		S.rr -= r_ - p_/(cr*cr);
	}

	if (ur<cr)
	{
		tmp = xl - (ur-cr)*dt;
		p_ = pr - get_CON_aveL(geom, xl, tmp, xr, p_par);
		u_ = ur - get_PRM_aveL(geom, xl, tmp, xr, u_par);

		B_ =  (u_*rr*cr - p_)/2.0;

		S.rr -= -B_/(cr*cr);
		S.pr -= -B_;
		S.ur -=  B_/(rr*cr);
	}

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
