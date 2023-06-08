
__device__ double get_PLM_aveR(int geom, double x, double* par)
{
	double sL = par[0];
	double aM = par[1];
	double sR = par[2];
	double val;

	val = aM-sL + (sR+sL)*((1.0+x)/2.0);

	return val;
}

__device__ double get_PLM_aveL(int geom, double x, double* par)
{
	double sL = par[0];
	double aM = par[1];
	double sR = par[2];
	double val;

	val = aM-sL + (sR+sL)*(x/2.0);

	return val;
}
<<<<<<< HEAD

__device__ double get_PLM_slopeR(int geom, double x1, double x2, double* par)
{
	double sL = par[0];
	double sR = par[2];
	double a6 = -3.0*(sR-sL);
	x2 = 1.0 - x2;
	x1 = 1.0 - x1;

	return (sR+sL-a6) + a6*(x1+x2)/1.5;
}

__device__ double get_PLM_slopeL(int geom, double x1, double x2, double* par)
{
	double sL = par[0];
	double sR = par[2];
	double a6 = -3.0*(sR-sL);

	return (sR+sL+a6) - a6*(x2+x1)/1.5;
}
=======
>>>>>>> methods
