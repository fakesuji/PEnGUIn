/*
__device__ void set_L_state(int i, int geom, double* xa, double* dx, double* dv, double rad, 
                            double* r, double* p, double* u, double* v, double* w, double dt, double us, State &S)
{
	double r_par[4], p_par[4], u_par[4];
	get_CON_parameters(i, geom, xa, dx, dv, r, r_par);
	get_CON_parameters(i, geom, xa, dx, dv, p, p_par);
	get_PRM_parameters(i, geom, xa, dx, dv, u, u_par);

	double xl = xa[i];
	double xr = xa[i+1];
	//double cl = sqrt(gam*p[i]/r[i]);
	//double tmp = xr - fmax(cl, u[i]+cl)*dt;
	double cl;
	double tmp = fmax(u[i],0.0);
	tmp = xr - sqrt(gam*p[i]/r[i] + tmp*tmp)*dt;

	if (tmp<xl) 
	{
		//printf("Error: reconstruction out of bound at L_state step 1 (%f, %f, %f), %f, %f, %f, %f, %f, \n", xl, tmp, xr, u[i-1], u[i], u[i+1], cl, dt); 
		tmp = xl;
	}

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

		Bp = -( (ul-u_) + (pl-p_)/C - us)/2.0/C;
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

		Bm =  ( (ul-u_) - (pl-p_)/C - us)/2.0/C;
	}

	S.pl = pl + (Bp + Bm)*C*C;
	S.ul = ul + (Bp - Bm)*C;
	S.rl = 1.0/( 1.0/rl - (B0 + Bp + Bm) );

	S.pl = fmax(smallp*p[i],S.pl);
	S.rl = fmax(smallr*r[i],S.rl);

	return;
}

__device__ void set_R_state(int i, int geom, double* xa, double* dx, double* dv, double rad, 
                            double* r, double* p, double* u, double* v, double* w, double dt, double us, State &S)
{
	double r_par[4], p_par[4], u_par[4];
	get_CON_parameters(i, geom, xa, dx, dv, r, r_par);
	get_CON_parameters(i, geom, xa, dx, dv, p, p_par);
	get_PRM_parameters(i, geom, xa, dx, dv, u, u_par);

	double xl = xa[i];
	double xr = xa[i+1];
	//double cr = sqrt(gam*p[i]/r[i]);
	//double tmp = xl + fmax(cr, -u[i]+cr)*dt;
	double cr;
	double tmp = fmax(-u[i],0.0);
	tmp = xl + sqrt(gam*p[i]/r[i] + tmp*tmp)*dt;

	if (tmp>xr) 
	{
		//printf("Error: reconstruction out of bound at R_state step 1 %f, %f, %f\n", xl, tmp, xr); 
		tmp = xr;
	}

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

		Bp = -( (ur-u_) + (pr-p_)/C - us)/2.0/C;
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

		Bm =  ( (ur-u_) - (pr-p_)/C - us)/2.0/C;
	}

	S.pr = pr + (Bp + Bm)*C*C;
	S.ur = ur + (Bp - Bm)*C;
	S.rr = 1.0/( 1.0/rr - (B0 + Bp + Bm) );

	S.pr = fmax(smallp*p[i],S.pr);
	S.rr = fmax(smallr*r[i],S.rr);

	return;
}
*/
__device__ void set_state(int i, int geom, double* xa, double* dx, double* dv, double rad, 
                          double* r, double* p, double* u, double* v, double* w, double dt, double us, State &S)
{
	double rl_par[4], pl_par[4], ul_par[4];
	double rr_par[4], pr_par[4], ur_par[4];

	double rl, pl, ul, cl;
	double rr, pr, ur, cr;
	double xl, xr, tmp, C;
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
	C = cl*rl;

	Bp = 0.0;
	Bm = 0.0;
	B0 = 0.0;

	if (ul>-cl)
	{
		tmp = xr - (ul+cl)*dt;

		//r_ = get_CON_aveR(geom, xl, tmp, xr, rl_par);
		p_ = get_CON_aveR(geom, xl, tmp, xr, pl_par);
		u_ = get_PRM_aveR(geom, xl, tmp, xr, ul_par);

		Bp = -( (ul-u_) + (pl-p_)/C - us)/2.0;
	}

	if (ul>0.0)
	{
		tmp = xr - ul*dt;

		r_ = get_CON_aveR(geom, xl, tmp, xr, rl_par);
		p_ = get_CON_aveR(geom, xl, tmp, xr, pl_par);

		B0 = (pl-p_)/C/C + 1.0/rl - 1.0/r_;
	}

	if (ul>cl)
	{
		tmp = xr - (ul-cl)*dt;

		//r_ = get_CON_aveR(geom, xl, tmp, xr, rl_par);
		p_ = get_CON_aveR(geom, xl, tmp, xr, pl_par);
		u_ = get_PRM_aveR(geom, xl, tmp, xr, ul_par);

		Bm =  ( (ul-u_) - (pl-p_)/C - us)/2.0;
	}

	S.rl = 1.0/( 1.0/rl - (B0 + Bp/C + Bm/C) );
	S.pl = pl + (Bp + Bm)*C;
	S.ul = ul + Bp - Bm;

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
	C = cr*rr;

	Bp = 0.0;
	Bm = 0.0;
	B0 = 0.0;

	if (ur<-cr)
	{
		tmp = xl - (ur+cr)*dt;

		//r_ = get_CON_aveL(geom, xl, tmp, xr, rr_par);
		p_ = get_CON_aveL(geom, xl, tmp, xr, pr_par);
		u_ = get_PRM_aveL(geom, xl, tmp, xr, ur_par);

		Bp = -( (ur-u_) + (pr-p_)/C - us)/2.0;
	}

	if (ur<0.0)
	{
		tmp = xl - ur*dt;

		r_ = get_CON_aveL(geom, xl, tmp, xr, rr_par);
		p_ = get_CON_aveL(geom, xl, tmp, xr, pr_par);

		B0 = (pr-p_)/C/C + 1.0/rr - 1.0/r_;
	}

	if (ur<cr)
	{
		tmp = xl - (ur-cr)*dt;

		//r_ = get_CON_aveL(geom, xl, tmp, xr, rr_par);
		p_ = get_CON_aveL(geom, xl, tmp, xr, pr_par);
		u_ = get_PRM_aveL(geom, xl, tmp, xr, ur_par);

		Bm =  ( (ur-u_) - (pr-p_)/C - us)/2.0;
	}

	S.rr = 1.0/( 1.0/rr - (B0 + Bp/C + Bm/C) );
	S.pr = pr + (Bp + Bm)*C;
	S.ur = ur + Bp - Bm;

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
