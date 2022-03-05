//=================================================================================

__device__ double get_pm_rare(State &S)
{
	double cl, cr, plz, prz, z, pm;

	z = gammfac;

	cl = sqrt(gam*S.pl/S.rl);
	cr = sqrt(gam*S.pr/S.rr);

	plz = pow(S.pl,z);
	prz = pow(S.pr,z);

	pm = pow( (cl + cr - gam*z*(S.ur-S.ul))/(cl/plz + cr/prz) , 1.0/z);
	pm = fmax(pm, smallp);

	return pm;
}

//=================================================================================

__device__ double get_pm_simple(State &S)
{
	double cl, cr, r_, c_, pm;

	cl = sqrt(gam*S.pl/S.rl);
	cr = sqrt(gam*S.pr/S.rr);

	r_ = 0.5*(S.rr+S.rl);
	c_ = 0.5*(cr+cl);

	pm = 0.5*(S.pr+S.pl) - 0.5*(S.ur-S.ul)*r_*c_;

	return pm;
}

//=================================================================================

__device__ double get_pm_iterative(State &S)
{
	double en_l, en_r, cl, cr, uml, umr;
	double pm, pl, rl, ul, pr, rr, ur, t1, t2;
	bool converge = false;
	int counter = 0;

	rl   = S.rl;
	rr   = S.rr;

	pl   = S.pl;
	pr   = S.pr;

	ul   = S.ul;
	ur   = S.ur;

	cl = sqrt(gam*pl*rl);
	cr = sqrt(gam*pr*rr);

	pm = (pl*cr + pr*cl - cl*cr*(ur-ul) ) / (cl + cr);

	while (!converge)
	{
		if (counter>16) return get_pm_rare(S); 

		t1   = (gampfac*pm + gammfac*pl) / pl;
		t2   = (gampfac*pm + gammfac*pr) / pr;
		t1   = cl * sqrt(t1);
		t2   = cr * sqrt(t2);
		uml  = ul - (pm - pl) / t1;
		umr  = ur + (pm - pr) / t2;

		en_l = 4.0 * t1 * t1 / rl;
		en_r = 4.0 * t2 * t2 / rr;
		en_l = -en_l * t1 / (en_l - gamp*(pm - pl));
		en_r =  en_r * t2 / (en_r - gamp*(pm - pr));

		t1   = (umr - uml)*(en_r * en_l) / (en_r - en_l);
		pm  += t1;

		converge = (fabs(t1)/pm<1.0e-12);
		counter++;

		if (signbit(pm)!=0)
		{
			converge = false;
			counter  = 100;
		}
  	}
	return pm;
}

//=================================================================================

__device__ double get_pm_iso(State &S)
{
	double rl, rr, du, pm;

	rl   = sqrt(S.rl);
	rr   = sqrt(S.rr);
	du   = S.ur-S.ul;

	double a, b, c;

	a = rl + rr;
	b = du*rl*rr/2.0;
	c = S.pr*rl + S.pl*rr;

	pm = (sqrt(b*b + a*c) - b) / a;

	return pm*pm;
}
//=================================================================================

__device__ void get_pm_um_rare(State &S, double &pm, double &um)
{
	double cl, cr, plz, prz, z;

	z = gammfac;

	cl = sqrt(gam*S.pl/S.rl);
	cr = sqrt(gam*S.pr/S.rr);

	plz = pow(S.pl,z);
	prz = pow(S.pr,z);

	pm = pow( (cl + cr - gam*z*(S.ur-S.ul))/(cl/plz + cr/prz) , 1.0/z);
	um = (plz*S.ul/cl + prz*S.ur/cr - (prz-plz)/z/gam) / (plz/cl + prz/cr);
  
	pm = fmax(pm, 0.0);

	return;
}

//=================================================================================

__device__ void get_pm_um_iso(State &S, double &pm, double &um)
{
	double rl, rr, du, t1, t2;

	rl   = sqrt(S.rl);
	rr   = sqrt(S.rr);
	du   = S.ur-S.ul;

	double a, b, c;

	a = rl + rr;
	b = du*rl*rr/2.0;
	c = S.pr*rl + S.pl*rr;

	pm = (sqrt(b*b + a*c) - b) / a;

	t1  = rl * pm;
	t2  = rr * pm;

	pm = pm*pm;

	um = 0.5 * (S.ur + S.ul + (pm-S.pr)/t2 - (pm-S.pl)/t1);

	return;
}

//=================================================================================

__device__ void get_pm_um_iterative(State &S, double &pm, double &um)
{
	double en_l, en_r, cl, cr, uml, umr;
	double pl, rl, ul, pr, rr, ur, t1, t2;
	bool converge = false;
	int counter = 0;

	rl = S.rl;
	rr = S.rr;

	pl = S.pl;
	pr = S.pr;

	ul = S.ul;
	ur = S.ur;

	cl = sqrt(gam*pl*rl);
	cr = sqrt(gam*pr*rr);

	pm = (pl*cr + pr*cl - cl*cr*(ur-ul) ) / (cl + cr);

	while (!converge)
	{
		if (counter>16)
		{
			get_pm_um_rare(S, pm, um);
			return; 
		}

		t1   = (gampfac*pm + gammfac*pl) / pl;
		t2   = (gampfac*pm + gammfac*pr) / pr;
		t1   = cl * sqrt(t1);
		t2   = cr * sqrt(t2);
		uml  = ul - (pm - pl) / t1;
		umr  = ur + (pm - pr) / t2;

		en_l = 4.0 * t1 * t1 / rl;
		en_r = 4.0 * t2 * t2 / rr;
		en_l = -en_l * t1 / (en_l - gamp*(pm - pl));
		en_r =  en_r * t2 / (en_r - gamp*(pm - pr));

		t1   = (umr - uml)*(en_r * en_l) / (en_r - en_l);
		pm  += t1;

		converge = (fabs(t1)/pm<1.0e-12);
		counter++;
    		if (signbit(pm)!=0)
		{
			converge = false;
			counter  = 100;
		}
	}
 	um = 0.5*(uml+umr);
	return;
}

//=================================================================================

__device__ void get_lr_speeds_direct(State S, double &sl, double &sr)
{
	double cl, cr, r_, c_, pm;

	cl = sqrt(gam*S.pl/S.rl);
	cr = sqrt(gam*S.pr/S.rr);

	r_ = 0.5*(S.rr+S.rl);
	c_ = 0.5*(cr+cl);

	pm = 0.5*(S.pr+S.pl) - 0.5*(S.ur-S.ul)*r_*c_;

	if (pm<=S.pl) sl = S.ul - cl;
	else          sl = S.ul - cl*sqrt(1.0 + ((gam+1.0)/(2.0*gam))*(pm/S.pl-1.0));

	if (pm<=S.pr) sr = S.ur + cr;
	else          sr = S.ur + cr*sqrt(1.0 + ((gam+1.0)/(2.0*gam))*(pm/S.pr-1.0));

	return;
}


//=================================================================================

__device__ void get_lr_speeds_Davis(State S, double &sl, double &sr)
{
	double hl, hr;

	hl = gam*S.pl/S.rl/gamm + 0.5*S.ul*S.ul;
	hr = gam*S.pr/S.rr/gamm + 0.5*S.ur*S.ur;

	double srl = sqrt(S.rl);
	double srr = sqrt(S.rr);

	double u_ = (srl*S.ul + srr*S.ur)/(srl+srr);
	double h_ = (srl*  hl + srr*  hr)/(srl+srr);

	double a_ = sqrt(gamm*(h_ - 0.5*u_*u_));

	sl = u_-a_;
	sr = u_+a_;

	return;
}

//=================================================================================

__device__ void get_lr_speeds_Einfeldt(State S, double &sl, double &sr)
{
	double cl, cr;

	cl = sqrt(gam*S.pl/S.rl);
	cr = sqrt(gam*S.pr/S.rr);

	double srl = sqrt(S.rl);
	double srr = sqrt(S.rr);

	double u_ = (srl*S.ul + srr*S.ur)/(srl+srr);

	double n2 = 0.5*srl*srr/((srl+srr)*(srl+srr));
	double c_ = sqrt( (cl*cl*srl + cr*cr*srr)/(srl+srr) + n2*(S.ur-S.ul)*(S.ur-S.ul) );

	sl = u_-c_;
	sr = u_+c_;

	return;
}

__device__ double get_sm(State S, double sl, double sr)
{
  return (S.pr-S.pl + S.rl*S.ul*(sl-S.ul) - S.rr*S.ur*(sr-S.ur))/(S.rl*(sl-S.ul) - S.rr*(sr-S.ur));
}

//=================================================================================

__device__ void wave_speeds(State S, double &pm, double &sl, double &sm, double &sr)
{
	#if EOS_flag>0
	get_pm_um_iterative(S, pm, sm);
	get_lr_speeds_Davis(S, sl, sr);
	//pm = get_pm_simple(S);
	//sm = get_sm(S,sl,sr);
	#else
	get_pm_um_iso(S, pm, sm);	
	get_lr_speeds_Einfeldt(S, sl, sr);
	#endif
	return;
}

//=================================================================================
