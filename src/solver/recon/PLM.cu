
__device__ double get_PLM_aveR(int geom, double x, double* par)
{
	double sR = par[2];
	double val;

	val = par[1] + sR*x;

	return val;
}

__device__ double get_PLM_aveL(int geom, double x, double* par)
{
	double sL = par[0];
	double val;

	val = par[1] - sL*(1.0-x);

	return val;
}
