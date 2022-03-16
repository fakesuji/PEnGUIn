
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
                         double u, double v, double w, body *planet, double dt)
{
	double gx, fx;

	#if geomx==0
		#ifdef shear_box
		gx = -rad;
		fx = 2.0*v;
		#else
		gx = 0.0;
		fx = 0.0;
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

		fx = v*v/rad;

	#elif geomx==2
		gx = -1.0/rad/rad;
		double gx_tmp, gy_tmp, gz_tmp;
		for (int m=0; m<n_planet; m++)
		{
			star_planet_grav(rad, azi, pol, planet[m], 0.0, gx_tmp, gy_tmp, gz_tmp);
			gx += gx_tmp;
		}

		fx = v*v/rad + w*w/rad;
	#endif
	return gx+fx;
}

//====================================================================

__device__ double get_fy(double rad, double azi, double pol,
                         double u, double v, double w, body *planet, double dt)
{
	double gy, fy;

	#if geomx==0
		#ifdef shear_box
		gy = 0.0;
		fy = 0.0;
		#else
		gy = 0.0;
		fy = 0.0;
		#endif

	#elif geomx==1
		gy = 0.0;
		double gx_tmp, gy_tmp, gz_tmp;
		for (int m=0; m<n_planet; m++)
		{
			star_planet_grav_cyl(rad, azi, pol, planet[m], 0.0, gx_tmp, gy_tmp, gz_tmp);
			gy += gy_tmp;
		}

		fy = 0.0;

	#elif geomx==2
		gy = 0.0;
		double gx_tmp, gy_tmp, gz_tmp;
		for (int m=0; m<n_planet; m++)
		{
			star_planet_grav(rad, azi, pol, planet[m], 0.0, gx_tmp, gy_tmp, gz_tmp);
			gy += gy_tmp;
		}

		fy = 0.0;
	#endif

	return gy+fy;
}

//====================================================================

__device__ double get_fz(double rad, double azi, double pol,
                         double u, double v, double w, body *planet, double dt)
{
	double gz, fz;

	#if geomx==0
		#ifdef shear_box
		gz = 0.0;
		fz = 0.0;
		#else
		gz = 0.0;
		fz = 0.0;
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
		
		fz = 0.0;

	#elif geomx==2
		gz = 0.0;
		double gx_tmp, gy_tmp, gz_tmp;
		for (int m=0; m<n_planet; m++)
		{
			star_planet_grav(rad, azi, pol, planet[m], 0.0, gx_tmp, gy_tmp, gz_tmp);
			gz += gz_tmp;
		}

		fz = v*v/rad/tan(pol);
		if (pol < zmin) {gz *= -1.0; fz *=-1.0;}
	#endif

	return gz+fz;
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
