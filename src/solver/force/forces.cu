__device__ double get_slope(int geom, double x0, double x1, double x2, double x3, double a0, double a1, double a2, double a3, double a4)
{
	//double d1 = x1-x0;
	double d2 = x2-x1;
	//double d3 = x3-x2;

	//return (a3-a1)/(0.5*d1+0.5*d3+d2);
	return (-(1.0/12.0)*a4+(2.0/3.0)*a3-(2.0/3.0)*a1+(1.0/12.0)*a0)/d2;
}


//====================================================================

__device__ void star_planet_grav(double rad, double azi, double pol, body p, double dt, double &fx, double &fy, double &fz)
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

	double Rp = sqrt(rad*rad + plx*plx - 2.0*plx*rad*cosfac*sinpol);
	double rs_fac = fmax(4.0-3.0*Rp/p.rs, 1.0)/fmax(p.rs*p.rs*p.rs, Rp*Rp*Rp);

	fx = -plm*(rad-plx*cosfac*sinpol)*rs_fac - plm*sinpol*cosfac/plx/plx;
	fy = -plm*plx*sinfac*rs_fac              + plm*sinpol*sinfac/plx/plx;
	fz =  plm*plx*cosfac*cospol*rs_fac       - plm*cospol*cosfac/(sinpol*plx*plx);
	return;
}

__device__ void star_planet_grav_cyl(double rad, double azi, double z, body p, double dt, double &fx, double &fy, double &fz)
{
	double cosfac, sinfac;
	sincos(azi-(p.y+dt*p.vy), &sinfac, &cosfac);

	double plm = p.m;
	double plx = p.x;

	double Rp = sqrt(rad*rad + plx*plx - 2.0*plx*rad*cosfac + z*z + p.rs*p.rs);

	fx = -plm*(rad-plx*cosfac)/(Rp*Rp*Rp) - plm*cosfac/plx/plx;
	fy = -plm*plx*sinfac/(Rp*Rp*Rp)       + plm*sinfac/plx/plx;
	fz = -plm*z/(Rp*Rp*Rp);
	return;
}

//====================================================================

__device__ double get_fx(double rad, double azi, double pol,
                         double u, double v, double w)
{
	double fx;

	#if geomx==0
		#ifdef shear_box
		fx = 2.0*v;
		#else
		fx = 0.0;
		#endif

	#elif geomx==1
		fx = v*v/rad;

	#elif geomx==2
		fx = v*v/rad + w*w/rad;
	#endif
	return fx;
}

//====================================================================

__device__ double get_fy(double rad, double azi, double pol,
                         double u, double v, double w)
{
	double fy;

	#if geomx==0
		#ifdef shear_box
		fy = 0.0;
		#else
		fy = 0.0;
		#endif

	#elif geomx==1
		fy = 0.0;

	#elif geomx==2
		fy = 0.0;
	#endif

	return fy;
}

//====================================================================

__device__ double get_fz(double rad, double azi, double pol,
                         double u, double v, double w)
{
	double fz;

	#if geomx==0
		#ifdef shear_box
		fz = 0.0;
		#else
		fz = 0.0;
		#endif

	#elif geomx==1
		fz = 0.0;

	#elif geomx==2
		fz = v*v/rad/tan(pol);
		if (pol < zmin) {fz *=-1.0;}
	#endif

	return fz;
}

//====================================================================

__device__ double get_gx(double rad, double azi, double pol, body *planet)
{
	double gx;

	#if geomx==0
		#ifdef shear_box
		gx = -rad;
		#else
		gx = 0.0;
		#endif

	#elif geomx==1
		double dis = sqrt(rad*rad+pol*pol);
		gx = -rad/dis/dis/dis;
		double gx_tmp, gy_tmp, gz_tmp;
		for (int m=0; m<n_planet; m++)
		{
			star_planet_grav_cyl(rad, azi, pol, planet[m], 0.0, gx_tmp, gy_tmp, gz_tmp);
			gx += gx_tmp;
		}

	#elif geomx==2
		gx = -1.0/rad/rad;
		double gx_tmp, gy_tmp, gz_tmp;
		for (int m=0; m<n_planet; m++)
		{
			star_planet_grav(rad, azi, pol, planet[m], 0.0, gx_tmp, gy_tmp, gz_tmp);
			gx += gx_tmp;
		}
	#endif
	return gx;
}

//====================================================================

__device__ double get_gy(double rad, double azi, double pol, body *planet)
{
	double gy;

	#if geomx==0
		#ifdef shear_box
		gy = 0.0;
		#else
		gy = 0.0;
		#endif

	#elif geomx==1
		gy = 0.0;
		double gx_tmp, gy_tmp, gz_tmp;
		for (int m=0; m<n_planet; m++)
		{
			star_planet_grav_cyl(rad, azi, pol, planet[m], 0.0, gx_tmp, gy_tmp, gz_tmp);
			gy += gy_tmp;
		}

	#elif geomx==2
		gy = 0.0;
		double gx_tmp, gy_tmp, gz_tmp;
		for (int m=0; m<n_planet; m++)
		{
			star_planet_grav(rad, azi, pol, planet[m], 0.0, gx_tmp, gy_tmp, gz_tmp);
			gy += gy_tmp;
		}
	#endif

	return gy;
}

//====================================================================

__device__ double get_gz(double rad, double azi, double pol, body *planet)
{
	double gz;

	#if geomx==0
		#ifdef shear_box
		gz = 0.0;
		#else
		gz = 0.0;
		#endif

	#elif geomx==1
		double dis = sqrt(rad*rad+pol*pol);
		gz = -pol/dis/dis/dis;
		double gx_tmp, gy_tmp, gz_tmp;
		for (int m=0; m<n_planet; m++)
		{
			star_planet_grav_cyl(rad, azi, pol, planet[m], 0.0, gx_tmp, gy_tmp, gz_tmp);
			gz += gz_tmp;
		}

	#elif geomx==2
		gz = 0.0;
		double gx_tmp, gy_tmp, gz_tmp;
		for (int m=0; m<n_planet; m++)
		{
			star_planet_grav(rad, azi, pol, planet[m], 0.0, gx_tmp, gy_tmp, gz_tmp);
			gz += gz_tmp;
		}
		if (pol < zmin) {gz *= -1.0;}
	#endif

	return gz;
}

//====================================================================

__device__ void get_grav(double rad, double azi, double pol,
                         body *planet, double dt,
                         double &gx, double &gy, double &gz)
{
	#if geomx==0
		#ifdef shear_box
		gx = -rad;
		gy = 0.0;
		gz = 0.0;
		#else
		gx = 0.0;
		gy = 0.0;
		gz = 0.0;
		#endif

	#elif geomx==1
		double dis = sqrt(rad*rad+pol*pol);
		gx = -rad/dis/dis/dis;
		gy = 0.0;
		gz = -pol/dis/dis/dis;

		double gx_tmp, gy_tmp, gz_tmp;
		for (int m=0; m<n_planet; m++)
		{
			star_planet_grav_cyl(rad, azi, pol, planet[m], 0.0, gx_tmp, gy_tmp, gz_tmp);
			gx += gx_tmp;
			gy += gy_tmp;
			gz += gz_tmp;
		}

	#elif geomx==2
		gx = -1.0/rad/rad;
		gy = 0.0;
		gz = 0.0;

		double gx_tmp, gy_tmp, gz_tmp;
		for (int m=0; m<n_planet; m++)
		{
			star_planet_grav(rad, azi, pol, planet[m], 0.0, gx_tmp, gy_tmp, gz_tmp);
			gx += gx_tmp;
			gy += gy_tmp;
			gz += gz_tmp;
		}
		if (pol < zmin) {gz *= -1.0;}
	#endif
	return;
}

//====================================================================

__device__ void get_fict(double rad, double azi, double pol,
                         double u, double v, double w, 
                         double &fx, double &fy, double &fz)
{
	#if geomx==0
		#ifdef shear_box
		fx = 2.0*v;
		fy = 0.0;
		fz = 0.0;
		#else
		fx = 0.0;
		fy = 0.0;
		fz = 0.0;
		#endif

	#elif geomx==1
		fx = v*v/rad;
		fy = 0.0;
		fz = 0.0;

	#elif geomx==2
		fx = v*v/rad + w*w/rad;
		fy = 0.0;
		fz = v*v/rad/tan(pol);
		if (pol < zmin) {fz *= -1.0;}
	#endif
	return;
}
