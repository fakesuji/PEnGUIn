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
//	label += "_test";

	printf("label %s assigned. \n\n", label.c_str());
 
	return label;
}

void write_grid_val(ofstream &ofile, Grid* G)
{
	double tmp;

	for (int k=zpad; k<zres+zpad; k++)
	for (int j=ypad; j<yres+ypad; j++)
	for (int n=0; n<ndev; n++)
	for (int i=xpad; i<G[n].xres+xpad; i++)
	{
		tmp = G[n].get_r(i,j,k);
		ofile.write((char*)&tmp, sizeof(double));		
	}

	for (int k=zpad; k<zres+zpad; k++)
	for (int j=ypad; j<yres+ypad; j++)
	for (int n=0; n<ndev; n++)
	for (int i=xpad; i<G[n].xres+xpad; i++)
	{
		tmp = G[n].get_p(i,j,k);
		ofile.write((char*)&tmp, sizeof(double));		
	}

	for (int k=zpad; k<zres+zpad; k++)
	for (int j=ypad; j<yres+ypad; j++)
	for (int n=0; n<ndev; n++)
	for (int i=xpad; i<G[n].xres+xpad; i++)
	{
		tmp = G[n].get_u(i,j,k);
		ofile.write((char*)&tmp, sizeof(double));		
	}

	#if ndim>1
	for (int k=zpad; k<zres+zpad; k++)
	for (int j=ypad; j<yres+ypad; j++)
	for (int n=0; n<ndev; n++)
	for (int i=xpad; i<G[n].xres+xpad; i++)
	{
		tmp = G[n].get_v(i,j,k);
		ofile.write((char*)&tmp, sizeof(double));		
	}
	#endif

	#if ndim>2
	for (int k=zpad; k<zres+zpad; k++)
	for (int j=ypad; j<yres+ypad; j++)
	for (int n=0; n<ndev; n++)
	for (int i=xpad; i<G[n].xres+xpad; i++)
	{
		tmp = G[n].get_w(i,j,k);
		ofile.write((char*)&tmp, sizeof(double));		
	}
	#endif

	#ifdef dust_flag

	for (int k=zpad; k<zres+zpad; k++)
	for (int j=ypad; j<yres+ypad; j++)
	for (int n=0; n<ndev; n++)
	for (int i=xpad; i<G[n].xres+xpad; i++)
	{
		tmp = G[n].get_r_dust(i,j,k);
		ofile.write((char*)&tmp, sizeof(double));		
	}

	for (int k=zpad; k<zres+zpad; k++)
	for (int j=ypad; j<yres+ypad; j++)
	for (int n=0; n<ndev; n++)
	for (int i=xpad; i<G[n].xres+xpad; i++)
	{
		tmp = G[n].get_u_dust(i,j,k);
		ofile.write((char*)&tmp, sizeof(double));		
	}

	#if ndim>1
	for (int k=zpad; k<zres+zpad; k++)
	for (int j=ypad; j<yres+ypad; j++)
	for (int n=0; n<ndev; n++)
	for (int i=xpad; i<G[n].xres+xpad; i++)
	{
		tmp = G[n].get_v_dust(i,j,k);
		ofile.write((char*)&tmp, sizeof(double));		
	}
	#endif

	#if ndim>2
	for (int k=zpad; k<zres+zpad; k++)
	for (int j=ypad; j<yres+ypad; j++)
	for (int n=0; n<ndev; n++)
	for (int i=xpad; i<G[n].xres+xpad; i++)
	{
		tmp = G[n].get_w_dust(i,j,k);
		ofile.write((char*)&tmp, sizeof(double));		
	}
	#endif

	#endif

	return;
}

void write_check_point(ofstream &ofile, double simtime, Grid* glo)
{
	ofile.write((char*)&simtime, sizeof(double));

	for (int i=xpad; i<xres+xpad+1; i++) ofile.write((char*)&glo[0].xa[i], sizeof(double));
	if (ndim>1) 
	for (int i=ypad; i<yres+ypad+1; i++) ofile.write((char*)&glo[0].ya[i], sizeof(double));
	if (ndim>2) 
	for (int i=zpad; i<zres+zpad+1; i++) ofile.write((char*)&glo[0].za[i], sizeof(double));

	//////////////////////////////////////////////

	write_grid_val(ofile, glo);

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
