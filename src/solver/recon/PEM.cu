__device__ double get_g_apprx(double x, double k)
{
	double y = x-1.0;
	double tmp = 1.0;
	double val = tmp;

	tmp *= y*(k+1.0)/2.0;
	val += tmp;
	
	tmp *= y*(k-1.0)/3.0;
	val += tmp;

	tmp *= y*(k-2.0)/4.0;
	val += tmp;

	return val;
}

//=======================================================================================

__device__ double get_PEM_0(double x, double s0, double aM, double eta, double lx)
{
	double val;
	val = aM + s0*(1.0-__expf(eta*lx));

	return val;
}

__device__ double get_PEM_1(double x, double aM, double s1, double eta, double lx)
{
	double val;
	double y=1.0-x;
	if (eta*y>1e-4) val = aM + s1*lim01(x*(1.0-__expf(eta*lx))/(y*eta));
	else            val = aM + s1*get_g_apprx(x, eta);

	return val;
}

//=======================================================================================

__device__ double get_PEM_aveR(int geom, double x, double* par, double lx, double ly)
{
	double sL = par[0];
	double aM = par[1];
	double sR = par[2];
	double val;

	if      (sR==0.0) val = aM;
	else if (sR/sL>=1.0) 
	{
		if      (x>=1.0) val = aM + sR;
		else if (x<=0.0) val = aM;
		else             val = get_PEM_1(x, aM, sR, sR/sL, lx);
	}
	else
	{
		x = 1.0-x;
		if      (x>=1.0) val = aM;
		else if (x<=0.0) val = aM + sR;
		else             val = get_PEM_0(x, sR, aM, sL/sR, ly);
	}

	return val;
}

__device__ double get_PEM_aveL(int geom, double x, double* par, double lx, double ly)
{
	double sL = par[0];
	double aM = par[1];
	double sR = par[2];
	double val;

	if      (sL==0.0) val = aM;
	else if (sR/sL>=1.0) 
	{
		if      (x>=1.0) val = aM;
		else if (x<=0.0) val = aM - sL;
		else             val = get_PEM_0(x, -sL, aM, sR/sL, lx);
	}
	else
	{
		x = 1.0-x;
		if      (x>=1.0) val = aM - sL;
		else if (x<=0.0) val = aM;
		else             val = get_PEM_1(x, aM, -sL, sL/sR, ly);
	}

	return val;
}
