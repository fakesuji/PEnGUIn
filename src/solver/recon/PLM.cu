
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
