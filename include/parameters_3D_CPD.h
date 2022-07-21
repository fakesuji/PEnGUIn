#ifndef PARAMETERS_H
#define PARAMETERS_H

#define dump_flag
#define kill_flag 1
//#define OrbAdv_flag
//#define visc_flag
//#define dust_flag

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
const double smallp = 1.0e-14;
const double smallr = 1.0e-14;

//=======================================================================
// Geometric parameters
//=======================================================================

#define recon_flag 2

const int std_thd = 1024;

#define ndim 3

#if recon_flag==2
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

const int xres = 192;
const int yres = 672;
const int zres = 64;

const int xarr = xres + 2*xpad;
const int yarr = yres + 2*ypad;
const int zarr = zres + 2*zpad;

const int x_xdiv = 16;
const int x_ydiv = 8;
const int x_zdiv = 1;

const int x_xthd = x_xdiv + 2*xpad;

const int y_xdiv = 8;
const int y_ydiv = 16;
const int y_zdiv = 1;

const int y_ythd = y_ydiv + 2*ypad;

const int z_xdiv = 8;
const int z_ydiv = 1;
const int z_zdiv = 16;

const int z_zthd = z_zdiv + 2*zpad;

//=======================================================================
// Geometric parameters
//=======================================================================

const int ndev = 1;

//=======================================================================
// Temporal parameters
//=======================================================================

const double frame_omega = 1.0;

const double sav_interval = 1.0*twopi;
const double sta_time = 0.0;
const double end_time = 20.0*twopi;

const int prt_interval = 100;
const int max_step = 1000000000;

//=======================================================================
// Hydro parameters
//=======================================================================

#define EOS_flag 0
#define internal_e_flag 0

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

const double p_beta = 7.0/7.0;                             // temperature ~ r^-p_beta
#if ndim==3
const double p_alpha = 1.5 - 0.5*p_beta + 1.5;          // midplane density ~ r^-p_alpha (isothermal limit)
#else 
const double p_alpha = 1.5;                             // surface density ~ r^-p_alpha
#endif 
const double ss_alpha = 0.004;                            // alpha-viscosity
const double sc_h = 0.035;                              // scale height at r=1, normalized to that at r = 100

#if ndim==3
const double Sigma_0 = 0.1*MMSN_1AU/(sqrt_hpi*sc_h);    // midplane density at r=1 in units of M_solar/AU^3
#else 
const double Sigma_0 = (MMSN_1AU/0.76)*0.05;              // density at r=1 in units of M_solar/AU^2
#endif

const double kill_width = 2.0;                      // in units of sc_h
const double beta_cool = 10.0;                          // in units of dynamical time for beta cooling

//=======================================================================
// planet parameters
//=======================================================================

const int n_planet = 1;
const double planet_mass = 3.0*0.00042875;
const double planet_radius = 1.0;
const double planet_ecc = 0.0;

const double ramp_time = twopi*5.0;

//=======================================================================
// Grid parameters
//=======================================================================

const double xmin = 0.65;
const double xmax = 1.35;
const double ymin = 0.0;
const double ymax = twopi;
const double zmin = hpi-4.0*sc_h;
const double zmax = hpi;

const double min_res = sc_h/32.0;
const double max_res = sc_h/3.0;
const double max_zres = sc_h/6.0;

#define geomx 2
#define geomy 4
#define geomz 5

const int gridx = 2;
const int gridy = 4;
const int gridz = 3;


#endif
