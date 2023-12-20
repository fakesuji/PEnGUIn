#ifndef PARAMETERS_H
#define PARAMETERS_H

#define dump_flag
#define kill_flag 1
//#define OrbAdv_flag
//#define visc_flag 1
//#define cool_flag 1
//#define twobd_flag 1
//#define silence_flag
//#define advec_flag
//#define ave_flag

//=======================================================================
// Save path
//=======================================================================

const char path_for_dump[]="/scratch/fung/";

//=======================================================================
// Constants
//=======================================================================

const double pi = 3.1415926535897932384626434;
const double hpi = 0.5*pi;
const double twopi = 2.0*pi;
const double sqrt_hpi = 1.2533141373155;
const double third = 1.0/3.0;
const double fourthd = 4.0/3.0;
const double EarthMass = 0.000003;
const double NeptuneMass = 0.00005;
const double SaturnMass = 0.0002857;
const double JupiterMass = 0.0009543;
const double MMSN_1AU = 0.00019126835;
const double smallp = 1.0e-10;
const double smallr = 1.0e-10;

//=======================================================================
// Geometric parameters
//=======================================================================

#define recon_flag 6

const int std_thd = 1024;

#define ndim 3

#if recon_flag>4
const int npad = 3;
#else
const int npad = 2;
#endif

const int xpad = npad;
#if ndim > 1
const int ypad = npad;
#else
const int ypad = 0;
#endif
#if ndim > 2
const int zpad = npad;
#else
const int zpad = 0;
#endif

const int xres = 384;
const int yres = 768;
const int zres = 216;

const int xarr = xres + 2*xpad;
const int yarr = yres + 2*ypad;
const int zarr = zres + 2*zpad;

//=======================================================================
// Domain division parameters
//=======================================================================

const int ndev = 1;

const int x_xdiv = 24;
const int x_ydiv = 8;
const int x_zdiv = 1;

const int x_xthd = x_xdiv + 2*xpad;

const int y_xdiv = 8;
const int y_ydiv = 24;
const int y_zdiv = 1;

const int y_ythd = y_ydiv + 2*ypad;

const int z_xdiv = 8;
const int z_ydiv = 1;
const int z_zdiv = 24;

const int z_zthd = z_zdiv + 2*zpad;

//=======================================================================
// Temporal parameters
//=======================================================================

const double frame_omega = 1.0;

const double sav_interval = 0.1*twopi/frame_omega;
const double end_time = 10.0*twopi/frame_omega;

const int prt_interval = 1000;
const int max_step = 1000000000;

//=======================================================================
// Hydro parameters
//=======================================================================

#define EOS_flag 2
#define internal_e_flag 1

#if EOS_flag == 0
const double gam = 1.0;
#else
const double gam = 1.4;
#endif
const double gamm = gam - 1.0;
const double gamp = gam + 1.0;
const double gammfac = gamm/gam/2.0;
const double gampfac = gamp/gam/2.0;
const double CFL = 0.3;

//=======================================================================
// Dust parameters
//=======================================================================

const double D_G_ratio = 0.01;
const double Stokes = 0.01;

//=======================================================================
// boundary parameters
//=======================================================================

#define init_flag 3

const int bound_lft = 0;
const int bound_rgh = 0;

const int bound_bak = 3;
const int bound_frn = 3;

const int bound_bom = 2;
const int bound_top = 2;

//=======================================================================
// Disk parameters
//=======================================================================

const double p_beta = 1.0;
#if ndim==3
const double p_alpha = 1.5 - 0.5*p_beta + 1.5;
#else 
const double p_alpha = 1.5;
#endif 
const double ss_alpha = 0.001;
const double sc_h = 0.05;

#if ndim==3
const double Sigma_0 = 0.1*MMSN_1AU/(sqrt_hpi*sc_h);
#else 
const double Sigma_0 = (MMSN_1AU/0.76)*0.05;
#endif

const double kill_width = 1.0;
const double beta_cool = 0.01;

//=======================================================================
// planet parameters
//=======================================================================

#if twobd_flag==1
const int n_planet = 2;
#else 
const int n_planet = 1;
#endif

const double planet_mass = 0.00003;
const double planet_radius = 1.0;
const double planet_ecc = 0.0;

const double rs_fac = 0.015;
const double ramp_time = 1.0*twopi/frame_omega;

//=======================================================================
// Grid parameters
//=======================================================================

const double xmin = 0.7;
const double xmax = 1.3;
const double ymin = 0.0;
const double ymax = twopi;
const double zmin = hpi-4.0*sc_h;
const double zmax = hpi;

#define geomx 2
#define geomy 4
#define geomz 5

const int gridx = 2;
const int gridy = 4;
const int gridz = 3;


#endif
