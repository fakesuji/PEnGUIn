//====================================================================

__device__ void o_star_planet_grav(double rad, double azi, double pol, body p, double dt, double &fx, double &fy, double &fz)
{
	double cosfac, sinfac;
	sincos(azi-(p.y+dt*p.vy), &sinfac, &cosfac);
	double cospol, sinpol;
	#if ndim==3
	sincos(pol, &sinpol, &cospol);
	#else
	cospol = 0.0;
	sinpol = 1.0;
	#endif

	double plm = p.m;
	double plx = p.x;
	
	double Rp, rs_fac, tmp;

	Rp = sqrt(rad*rad + plx*plx - 2.0*plx*rad*cosfac*sinpol);
	rs_fac = fmax(4.0-3.0*Rp/p.rs, 1.0)/fmax(p.rs*p.rs*p.rs, Rp*Rp*Rp);
	if (Rp<p.rs)
	{
		tmp = sin(hpi*Rp/p.rs);
		rs_fac *= tmp*tmp;
	}

	//double Rp = sqrt(rad*rad + plx*plx - 2.0*plx*rad*cosfac*sinpol + p.rs*p.rs);
	//double rs_fac = 1.0/Rp/Rp/Rp;

	fx = -plm*(rad-plx*cosfac*sinpol)*rs_fac - plm*sinpol*cosfac/plx/plx;
	fy = -plm*plx*sinfac*rs_fac              + plm*sinpol*sinfac/plx/plx;
	fz =  plm*plx*cosfac*cospol*rs_fac       - plm*cospol*cosfac/(sinpol*plx*plx);
	return;
}

__device__ void o_star_planet_grav_cyl(double rad, double azi, double z, body p, double dt, double &fx, double &fy, double &fz)
{
	double cosfac, sinfac;
	sincos(azi-(p.y+dt*p.vy), &sinfac, &cosfac);

	double plm = p.m;
	double plx = p.x;

	double Rp = sqrt(rad*rad + plx*plx - 2.0*plx*rad*cosfac + z*z + p.rs*p.rs);
	double rs_fac = 1.0/Rp/Rp/Rp;
	//double Rp = sqrt(rad*rad + plx*plx - 2.0*plx*rad*cosfac);
	//double rs_fac = fmax(4.0-3.0*Rp/p.rs, 1.0)/fmax(p.rs*p.rs*p.rs, Rp*Rp*Rp);

	fx = -plm*(rad-plx*cosfac)*rs_fac - plm*cosfac/plx/plx;
	fy = -plm*plx*sinfac*rs_fac       + plm*sinfac/plx/plx;
	fz = -plm*z*rs_fac;
	return;
}

//====================================================================

__device__ void o_star_planet_grav_inertial(double rad, double azi, double pol, body p, double dt, double &fx, double &fy, double &fz)
{
	double cosfac, sinfac;
	sincos(azi-(p.y+dt*p.vy), &sinfac, &cosfac);
	double cospol, sinpol;
	#if ndim==3
	sincos(pol, &sinpol, &cospol);
	#else
	cospol = 0.0;
	sinpol = 1.0;
	#endif

	double plm = p.m;
	double plx = p.x;

	//double Rp = sqrt(rad*rad + plx*plx - 2.0*plx*rad*cosfac*sinpol);
	//double rs_fac = fmax(4.0-3.0*Rp/p.rs, 1.0)/fmax(p.rs*p.rs*p.rs, Rp*Rp*Rp);

	double Rp = sqrt(rad*rad + plx*plx - 2.0*plx*rad*cosfac*sinpol + p.rs*p.rs);
	double rs_fac = 1.0/Rp/Rp/Rp;

	fx = -plm*(rad-plx*cosfac*sinpol)*rs_fac;
	fy = -plm*plx*sinfac*rs_fac;
	fz =  plm*plx*cosfac*cospol*rs_fac;
	return;
}

__device__ void o_star_planet_grav_cyl_inertial(double rad, double azi, double z, body p, double dt, double &fx, double &fy, double &fz)
{
	double cosfac, sinfac;
	sincos(azi-(p.y+dt*p.vy), &sinfac, &cosfac);

	double plm = p.m;
	double plx = p.x;

	double Rp = sqrt(rad*rad + plx*plx - 2.0*plx*rad*cosfac + z*z + p.rs*p.rs);
	double rs_fac = 1.0/Rp/Rp/Rp;
	//double Rp = sqrt(rad*rad + plx*plx - 2.0*plx*rad*cosfac);
	//double rs_fac = fmax(4.0-3.0*Rp/p.rs, 1.0)/fmax(p.rs*p.rs*p.rs, Rp*Rp*Rp);

	fx = -plm*(rad-plx*cosfac)*rs_fac;
	fy = -plm*plx*sinfac*rs_fac;
	fz = -plm*z*rs_fac;
	return;
}

//====================================================================

__device__ double output_gy(double rad, double azi, double pol, body planet)
{
	double gy = 0.0;
	double gx_tmp, gz_tmp;

	#if geomx==1
		#if twobd_flag == 1
		o_star_planet_grav_cyl_inertial(rad, azi, pol, planet, 0.0, gx_tmp, gy, gz_tmp);
		#else
		o_star_planet_grav_cyl(rad, azi, pol, planet, 0.0, gx_tmp, gy, gz_tmp);
		#endif

	#elif geomx==2
		#if twobd_flag == 1
		o_star_planet_grav_inertial(rad, azi, pol, planet, 0.0, gx_tmp, gy, gz_tmp);
		#else
		o_star_planet_grav(rad, azi, pol, planet, 0.0, gx_tmp, gy, gz_tmp);
		#endif
	#endif

	return gy;
}
