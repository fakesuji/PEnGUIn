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

__device__ double get_PEM_aveL(int geom, double x, double* par)
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
}


