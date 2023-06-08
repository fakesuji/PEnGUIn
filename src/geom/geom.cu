#include "parameters.h"

double get_dv_dr(int geom, double ra, double dr)
{
	if 	(geom==5)
	{
		if (dr<1.0e-4) return sin(ra+0.5*dr);
		else           return (cos(ra)-cos(ra+dr))/dr;
	}
	//else if (geom==4) return 1.0;	
	//else if (geom==3) return 1.0;
	else if (geom==2) return ra*ra + ra*dr + dr*dr/3.0;
	else if (geom==1) return ra + dr/2.0;
	else return 1.0;
}

double get_volume(int geom, double ra, double dr)
{
	return dr*get_dv_dr(geom,ra,dr);
}

void linspace(double* a, double start, double end, int len)
{
	for (int i=0; i<len; i++) a[i] = start + (end-start)*(double)i/(double)(len-1);
	return;
}

void logspace(double* a, double start, double end, int len)
{
	double tmp = log(end/start);
	for (int i=0; i<len; i++) a[i] = start*exp(tmp*(double)i/(double)(len-1));//pow(end/start,(double)i/(double)(len-1))*start;
	return;
}

void nonuspace(double* a, double start, double end, int len)
{
	int N = (len-1)/2;
	double L = (end-start)/2.0;

<<<<<<< HEAD
	double dx_min = min_res;
	double dx_max = max_res;
=======
	double dx_min = (L/(double)N)/2.0;
	double dx_max = (L/(double)N)*2.0;
>>>>>>> methods

	double k = log(1.0 - (dx_max-dx_min)/(L-(double)N*dx_min))/log(1.0-1.0/(double)N);

	for (int i=0; i<N; i++) a[i] = (end+start)/2.0 - (N-i)*dx_min - (L - N*dx_min)*pow((double)(N-i)/(double)N, k);
	for (int i=0; i<N+1; i++) a[N+i] = (end+start)/2.0 + i*dx_min + (L - N*dx_min)*pow((double)i/(double)N, k);

	return;
}

<<<<<<< HEAD
=======
void nonuspace_mix(double* a, double start, double end, int len)
{
	int N = xres/2;
	int N2 = (len-1-2*N)/2;
	double L = (xmax-xmin)/2.0;
	double L2 = (end-start-2.0*L)/2.0;

	double dx_min = (L/(double)N)/2.0;
	double dx_max = L2/(double)N2;

	double k = log(1.0 - (dx_max-dx_min)/(L-(double)N*dx_min))/log(1.0-1.0/(double)N);

	for (int i=0; i<N2; i++) a[i] = start + (double)i*dx_max;

	for (int i=0; i<N; i++) a[i+N2] = (end+start)/2.0 - (N-i)*dx_min - (L - N*dx_min)*pow((double)(N-i)/(double)N, k);

	for (int i=0; i<N; i++) a[i+N2+N] = (end+start)/2.0 + i*dx_min + (L - N*dx_min)*pow((double)i/(double)N, k);

	for (int i=0; i<N2+1; i++) a[i+N2+2*N] = start + 2.0*L + L2 + (double)i*dx_max;

	return;
}

>>>>>>> methods
void nonuspace_half(double* a, double start, double end, int len)
{
	int N = len-1;
	double L = end-start;

<<<<<<< HEAD
	double dx_min = min_res;
	double dx_max = max_res;
=======
	double dx_min = (L/(double)N)/2.0;
	double dx_max = (L/(double)N)*2.0;
>>>>>>> methods

	double k = log(1.0 - (dx_max-dx_min)/(L-(double)N*dx_min))/log(1.0-1.0/(double)N);

	for (int i=0; i<N+1; i++) a[i] = end - (N-i)*dx_min - (L - N*dx_min)*pow((double)(N-i)/(double)N, k);

	return;
}
