#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

#include "parameters.h"
#include "structs.h"
#include "util.h"

string create_label()
{
	string label;
	if      (ndim==1) label = int_to_string(xres);
	else if (ndim==2) label = int_to_string(xres)+"x"+int_to_string(yres);
	else if (ndim==3) label = int_to_string(xres)+"x"+int_to_string(yres)+"x"+int_to_string(zres);
	
	label += "_h"+int_to_string(1000.0*sc_h);

	if (n_planet>0)
	{
		if (planet_mass>=JupiterMass)
		{
			label += "_"+int_to_string(n_planet)+"p"+int_to_string(planet_mass/JupiterMass)+"J";
		}
		else
		{
			label += "_"+int_to_string(n_planet)+"p"+int_to_string(planet_mass/EarthMass)+"E";
		}
		label += "_e" +int_to_string(planet_ecc*1000.0);
	}

	#if visc_flag == 1
	label += "_a"+int_to_string(log10(ss_alpha)*10.0);
	#elif visc_flag == 2
	label += "_ah"+int_to_string(log10(ss_alpha)*10.0);
	#endif

	#ifdef cool_flag
	label += "_b"+int_to_string(log10(beta_cool)*10.0);
	#endif
	
	#ifdef OrbAdv_flag
	label += "_OA";
	#endif

	#ifdef dust_flag
	label += "_St"+int_to_string(log10(Stokes)*10.0);
	#endif

	#if recon_flag == 0
	label += "_VAN";
	#elif recon_flag == 1
	label += "_MOC";
	#elif recon_flag == 2
	label += "_PEM2";
	#elif recon_flag == 3
	label += "_PEM3";
	#elif recon_flag == 4
	label += "_PPM3";
	#elif recon_flag == 5
	label += "_PEM4";
	#elif recon_flag == 6
	label += "_PPM4";
	#endif

	#if init_flag == 6
	label += "_CFL"+int_to_string(CFL*100.0);
	#endif

	#ifdef rev_flag
	label += "_rev";
	#endif

	//label += "_"+int_to_string(ndev)+"dev";
	//label += "_video";

	printf("label %s assigned. \n\n", label.c_str());
 
	return label;
}


void restructure_data(double* data, Grid* G, int data_type)
{
	int ii, jj, kk, glo_idx,loc_idx;

	for (int k=zpad; k<zres+zpad; k++)
	for (int j=ypad; j<yres+ypad; j++)
	for (int n=0; n<ndev; n++)
	for (int i=xpad; i<G[n].xres+xpad; i++)
	{
		ii = i + G[n].xbgn - xpad;
		jj = j - ypad;
		kk = k - zpad;
		glo_idx = ii + xres*(jj + yres*kk);
		loc_idx = i + G[n].xarr*(j + G[n].yarr*kk);
		if      (data_type==0) data[glo_idx] = G[n].C[loc_idx].r;
		else if (data_type==1) data[glo_idx] = G[n].C[loc_idx].p;
		else if (data_type==2) data[glo_idx] = G[n].C[loc_idx].u;
		else if (data_type==3) data[glo_idx] = G[n].C[loc_idx].v;
		else if (data_type==4) data[glo_idx] = G[n].C[loc_idx].w;
		#ifdef dust_flag
		else if (data_type==5) data[glo_idx] = G[n].D[loc_idx].r;
		else if (data_type==6) data[glo_idx] = G[n].D[loc_idx].u;
		else if (data_type==7) data[glo_idx] = G[n].D[loc_idx].v;
		else if (data_type==8) data[glo_idx] = G[n].D[loc_idx].w;
		#endif
		//data[glo_idx] = 0.0;
	}
	return;
}

void write_check_point(ofstream &ofile, double simtime, Grid* G)
{
	size_t memsize= xres*yres*zres*sizeof(double);
	double* tmp = new double[xres*yres*zres];

	ofile.write((char*)&simtime, sizeof(double));

	ofile.write((char*)&G[0].xa[xpad], sizeof(double)*(xres+1));
	#if ndim>1
	ofile.write((char*)&G[0].ya[ypad], sizeof(double)*(yres+1));
	#endif
	#if ndim>2
	ofile.write((char*)&G[0].za[zpad], sizeof(double)*(zres+1));
	#endif

	restructure_data(tmp, G, 0);
	ofile.write((char*)&tmp[0],memsize);

	restructure_data(tmp, G, 1);
	ofile.write((char*)&tmp[0],memsize);

	restructure_data(tmp, G, 2);
	ofile.write((char*)&tmp[0],memsize);

	#if ndim>1
	restructure_data(tmp, G, 3);
	ofile.write((char*)&tmp[0],memsize);
	#endif

	#if ndim>2
	restructure_data(tmp, G, 4);
	ofile.write((char*)&tmp[0],memsize);
	#endif

	#ifdef dust_flag
	restructure_data(tmp, G, 5);
	ofile.write((char*)&tmp[0],memsize);

	restructure_data(tmp, G, 6);
	ofile.write((char*)&tmp[0],memsize);

	#if ndim>1
	restructure_data(tmp, G, 7);
	ofile.write((char*)&tmp[0],memsize);
	#endif

	#if ndim>1
	restructure_data(tmp, G, 8);
	ofile.write((char*)&tmp[0],memsize);
	#endif
	#endif

	delete[] tmp;

	return;
}

double load_grid(Grid* G, string fname)
{
	ifstream start_point;
	open_binary_file(start_point,fname);

	double start_time;
	start_point.read((char*)&start_time, sizeof(double));
	if (ndim>2) 
		start_point.seekg((1+xres+1+yres+1+zres+1)*sizeof(double), ios::beg);
	else if (ndim>1) 
		start_point.seekg((1+xres+1+yres+1)*sizeof(double), ios::beg);
	else
		start_point.seekg((1+xres+1)*sizeof(double), ios::beg);

	
	double tmp;
	for (int k=zpad; k<zres+zpad; k++)
	for (int j=ypad; j<yres+ypad; j++)
	for (int n=0; n<ndev; n++)
	for (int i=xpad; i<G[n].xres+xpad; i++)
	{
		start_point.read((char*)&tmp, sizeof(double));	
		G[n].write_r(i,j,k,tmp);	
	}

	for (int k=zpad; k<zres+zpad; k++)
	for (int j=ypad; j<yres+ypad; j++)
	for (int n=0; n<ndev; n++)
	for (int i=xpad; i<G[n].xres+xpad; i++)
	{
		start_point.read((char*)&tmp, sizeof(double));	
		G[n].write_p(i,j,k,tmp);	
	}

	for (int k=zpad; k<zres+zpad; k++)
	for (int j=ypad; j<yres+ypad; j++)
	for (int n=0; n<ndev; n++)
	for (int i=xpad; i<G[n].xres+xpad; i++)
	{
		start_point.read((char*)&tmp, sizeof(double));	
		G[n].write_u(i,j,k,tmp);	
	}

	if (ndim>1) 
	for (int k=zpad; k<zres+zpad; k++)
	for (int j=ypad; j<yres+ypad; j++)
	for (int n=0; n<ndev; n++)
	for (int i=xpad; i<G[n].xres+xpad; i++)
	{
		start_point.read((char*)&tmp, sizeof(double));	
		G[n].write_v(i,j,k,tmp);	
	}

	if (ndim>2) 
	for (int k=zpad; k<zres+zpad; k++)
	for (int j=ypad; j<yres+ypad; j++)
	for (int n=0; n<ndev; n++)
	for (int i=xpad; i<G[n].xres+xpad; i++)
	{
		start_point.read((char*)&tmp, sizeof(double));	
		G[n].write_w(i,j,k,tmp);	
	}

	close_output_file(start_point);

	return start_time;
}
