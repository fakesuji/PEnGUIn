//=======================================================================================

__device__ double get_PEM_0(double x, double S, double aB, double eta)
{
	double val;

	if (x==0.0) 
	{
		val = 0.0;
	}
	else
	{
	        val = exp(eta*__logf(x));
	}

	return aB + S*val;
}

__device__ double get_PEM_1(double x, double S, double aB, double eta)
{
	double tmp, val;

	if (x>=1.0)
	{
		val = 1.0;
	}
	if (eta*x>0.001)
	{
		val = (1.0-x)*(exp(eta*__logf(1.0-x))-1.0)/(eta*x)+1.0;
	}
	else
	{
		tmp = (eta+1.0)*x/2.0;
       		val = tmp;
        
        	tmp*=-(eta-1.0)*x/3.0;
        	val+= tmp;
        
        	tmp*=-(eta-2.0)*x/4.0;
        	val+= tmp;
        
        	tmp*=-(eta-3.0)*x/5.0;
        	val+= tmp;
	}

	return aB + S*val;
}

//=======================================================================================

<<<<<<< HEAD
__device__ double get_PEM_0(double x, double s0, double aM, double eta, double lx)
{
	double val;
	val = aM + s0*(1.0-__expf(eta*__logf(x)));

	return val;
}

__device__ double get_PEM_1(double x, double aM, double s1, double eta, double lx)
{
	double val;
	double y=1.0-x;
	if (eta*y>1e-4) val = aM + s1*lim01(x*(1.0-__expf(eta*__logf(x)))/(y*eta));
	else            val = aM + s1*get_g_apprx(x, eta);

	return val;
}

//=======================================================================================

__device__ double get_PEM_aveR(int geom, double x, double* par, double lx, double ly)
=======
__device__ double get_PEM_aveL(int geom, double x, double* par)
>>>>>>> methods
{
	double sL = par[0];
	double aM = par[1];
	double sR = par[2];

	if      (sL*sR==0.0) 
	{
		return aM;
	}
	else if (sR/sL>=1.0) 
	{
		if (sR/sL<=2.0) return get_PPM_aveL(geom, x, par);
		else            return get_PEM_0( x, sL, aM-sL, sR/sL);
	}
	else
	{
		if (sL/sR<=2.0) return get_PPM_aveL(geom, x, par);
		else            return get_PEM_1( x, sL, aM-sL, sL/sR);
	}
}

__device__ double get_PEM_aveR(int geom, double x, double* par)
{
	double sL = par[0];
	double aM = par[1];
	double sR = par[2];

	if      (sR*sL==0.0)
	{
		return aM;
	}
<<<<<<< HEAD
	else
	{
		x = 1.0-x;
		if      (x>=1.0) val = aM - sL;
		else if (x<=0.0) val = aM;
		else             val = get_PEM_1(x, aM, -sL, sL/sR, ly);
	}

	return val;
}

__device__ double get_PEM_slopeR(int geom, double x1, double x2, double* par)
{
	double sL = par[0];
	double sR = par[2];
/*
	if (sR/sL>=1.0) 
	{
		return (sR+sL)*(__powf(x2,sR/sL) - __powf(x1,sR/sL))/(x2-x1);
	}
	else
	{
		x2 = 1.0 - x2;
		x1 = 1.0 - x1;

		return (-sR-sL)*(__powf(x2,sL/sR) - __powf(x1,sL/sR))/(x2-x1);
	}
*/
	double slope = (4.0*sR-2.0*sL) - 2.0*(sR-sL)*(2.0-x2-x1);

	if (sR>=0.0) return fmax(slope,0.0);
	else         return fmin(slope,0.0);

	//return 2.0*(get_PEM_aveR(geom, x2, par, 0.0,0.0)-get_PEM_aveR(geom, x1, par, 0.0,0.0))/(x2-x1);
}

__device__ double get_PEM_slopeL(int geom, double x1, double x2, double* par)
{
	double sL = par[0];
	double sR = par[2];
/*
	if (sR/sL>=1.0) 
	{
		return (sR+sL)*(__powf(x2,sR/sL) - __powf(x1,sR/sL))/(x2-x1);
	}
	else
	{
		x2 = 1.0 - x2;
		x1 = 1.0 - x1;

		return (-sR-sL)*(__powf(x2,sL/sR) - __powf(x1,sL/sR))/(x2-x1);
	}
*/
	double slope = (4.0*sL-2.0*sR) + 2.0*(sR-sL)*(x2+x1);

	if (sR>=0.0) return fmax(slope,0.0);
	else         return fmin(slope,0.0);

	//return 2.0*(get_PEM_aveL(geom, x2, par, 0.0,0.0)-get_PEM_aveL(geom, x1, par, 0.0,0.0))/(x2-x1);
=======
	else if (sR/sL>1.0) 
	{
		if (sR/sL<=2.0) return get_PPM_aveR(geom, x, par);
		else            return get_PEM_1( 1.0-x, -sR, aM+sR, sR/sL);
	}
	else
	{
		if (sL/sR<=2.0) return get_PPM_aveR(geom, x, par);
		else            return get_PEM_0( 1.0-x, -sR, aM+sR, sL/sR);
	}
>>>>>>> methods
}


