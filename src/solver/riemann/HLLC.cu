
__device__ double get_flux(double FL, double F_L, double F_R, double FR, double sl, double sm, double sr)
{
	if (signbit(sl)==0)      return FL;
	else if (signbit(sr)!=0) return FR;
	else if (signbit(sm)==0) return F_L;
	else                     return F_R;
}

__device__ void HLLC_fluxes(State S, double pm, double sl, double sm, double sr, Cell &F, double &pres, double &uprs)
{
	double fl, f_l, f_r, fr;

	/////////////////////////////////////////////////////////

	fl  = S.ul*S.rl;
	f_l = sm*S.rl*(sl-S.ul)/(sl-sm);
	f_r = sm*S.rr*(sr-S.ur)/(sr-sm);
	fr  = S.ur*S.rr;
	F.r = get_flux(fl, f_l, f_r, fr, sl, sm, sr);

	F.u = get_flux(S.ul, S.ul, S.ur, S.ur, sl, sm, sr)*F.r;

	fl = S.vl;
	fr = S.vr;
	F.v = get_flux(fl, fl, fr, fr, sl, sm, sr)*F.r;

	fl = S.wl;
	fr = S.wr;
	F.w = get_flux(fl, fl, fr, fr, sl, sm, sr)*F.r;

	#if EOS_flag > 0
	fl = get_energy(S.rl, S.pl, S.ul, 0.0, 0.0);
	fr = get_energy(S.rr, S.pr, S.ur, 0.0, 0.0);
	F.p = get_flux(fl, fl, fr, fr, sl, sm, sr)*F.r;
	#else
	F.p = 0.0;
	#endif

	/////////////////////////////////////////////////////////

	fl  = S.pl;
	f_l = (sl*pm-sm*S.pl)/(sl-sm);
	f_r = (sr*pm-sm*S.pr)/(sr-sm);
	fr  = S.pr;
	pres = get_flux(fl, f_l, f_r, fr, sl, sm, sr);


	#if EOS_flag==2

	#if internal_e_flag==0
	fl  = S.ul*S.pl;
	f_l = sm*(sl*pm-S.ul*S.pl)/(sl-sm);
	f_r = sm*(sr*pm-S.ur*S.pr)/(sr-sm);
	fr  = S.ur*S.pr;
	#elif internal_e_flag==1
	fl  = S.ul;
	f_l = sm*(sl-S.ul)/(sl-sm);
	f_r = sm*(sr-S.ur)/(sr-sm);
	fr  = S.ur;
	#endif

	uprs = get_flux(fl, f_l, f_r, fr, sl, sm, sr);

	#else

	uprs = 0.0;

	#endif

	return;
}

__device__ void HLLE_fluxes(State S, double pm, double sl, double sm, double sr, Cell &F, double &pres, double &uprs)
{
	double fl, fr;

	/////////////////////////////////////////////////////////

	//sl = fmin(0.0,sl);
	//sr = fmax(0.0,sr);

	fl  = S.ul*S.rl;
	fr  = S.ur*S.rr;
	F.r = (sr*fl - sl*fr)/(sr-sl) + sr*sl*(S.rr-S.rl)/(sr-sl);

	fl  = S.ul*S.ul*S.rl;
	fr  = S.ur*S.ur*S.rr;
	F.u = (sr*fl - sl*fr)/(sr-sl) + sr*sl*(S.ur*S.rr-S.ul*S.rl)/(sr-sl);

	fl  = S.ul*S.vl*S.rl;
	fr  = S.ur*S.vr*S.rr;
	F.v = (sr*fl - sl*fr)/(sr-sl) + sr*sl*(S.vr-S.vl)/(sr-sl);

	fl  = S.ul*S.wl*S.rl;
	fr  = S.ur*S.wr*S.rr;
	F.w = (sr*fl - sl*fr)/(sr-sl) + sr*sl*(S.wr-S.wl)/(sr-sl);

	#if EOS_flag > 0
	double el = get_energy(S.rl, S.pl, S.ul, 0.0, 0.0);
	double er = get_energy(S.rr, S.pr, S.ur, 0.0, 0.0);

	fl  = S.ul*el*S.rl;
	fr  = S.ur*er*S.rr;
	F.p = (sr*fl - sl*fr)/(sr-sl) + sr*sl*(er*S.rr-el*S.rl)/(sr-sl);
	#endif

	/////////////////////////////////////////////////////////

	pres = (sr*S.pl-sl*S.pr)/(sr-sl);
	#if EOS_flag==2
	#if internal_e_flag==0
	uprs = (sr*S.ul*S.pl-sl*S.ur*S.pr)/(sr-sl);
	#elif internal_e_flag==1
	uprs = (sr*S.ul-sl*S.ur)/(sr-sl);
	#endif
	#else
	uprs = 0.0;
	#endif

	return;
}
