__host__ __device__ double get_dv_dr(int geom, double ra, double dr)
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

__host__ __device__ double get_volume(int geom, double ra, double dr)
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
