#ifndef PARAMETERS_H
#define PARAMETERS_H

#define mode_flag 0
#define dump_flag
//#define kill_flag 3
#define OrbAdv_flag
#define visc_flag
//#define advec_flag

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

const int std_thd = 1024;

#define ndim 2

const int npad = 2;

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

const int xres = 480;
const int yres = 1008;
const int zres = 1;

const int xarr = xres + 2*xpad;
const int yarr = yres + 2*ypad;
const int zarr = zres + 2*zpad;

const int x_xdiv = 48;
const int x_ydiv = 8;
const int x_zdiv = 1;

const int x_xthd = x_xdiv + 2*xpad;

const int y_xdiv = 16;
const int y_ydiv = 24;
const int y_zdiv = 1;

const int y_ythd = y_ydiv + 2*ypad;

const int z_xdiv = 16;
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

const double frame_omega = 42.5872131567;

const double sav_interval = 1.25*twopi/frame_omega;
const double end_time = 1000.0*twopi/frame_omega;

const int prt_interval = 1000;
const int max_step = 1000000000;

//=======================================================================
// Hydro parameters
//=======================================================================

#define recon_flag 1

#define EOS_flag 0
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
const double CFL = 0.2;

//=======================================================================
// boundary parameters
//=======================================================================

#define init_flag 2

const int bound_lft = 1;
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
const double p_alpha = 2.24 - 0.5*p_beta + 1.5;          // midplane density ~ r^-p_alpha (isothermal limit)
#else 
const double p_alpha = 2.24;                             // surface density ~ r^-p_alpha
#endif 
const double ss_alpha = 0.01;                            // alpha-viscosity
const double sc_h = 0.025;                              // scale height at r=1, normalized to that at r = 100

#if ndim==3
const double Sigma_0 = 0.1*MMSN_1AU/(sqrt_hpi*sc_h);    // midplane density at r=1 in units of M_solar/AU^3
#else 
const double Sigma_0 = (MMSN_1AU/0.76)*0.05;              // density at r=1 in units of M_solar/AU^2
#endif

const double kill_width = 4.0;                      // in units of sc_h
const double beta_cool = 10.0;                          // in units of dynamical time for beta cooling

//=======================================================================
// planet parameters
//=======================================================================

const int n_planet = 1;
const double planet_mass = 0.0123;
const double planet_radius = 0.082;
const double planet_ecc = 0.25;

const double ramp_time = 10.0*twopi/frame_omega;

//=======================================================================
// Grid parameters
//=======================================================================

const double xmin = 0.02;
const double xmax = 0.4;
const double ymin = 0.0;
const double ymax = twopi;
const double zmin = hpi-4.0*sc_h;
const double zmax = hpi;

#define geomx 1
#define geomy 3
#define geomz 0

const int gridx = 1;
const int gridy = 0;
const int gridz = 0;


#endif
