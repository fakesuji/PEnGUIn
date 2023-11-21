"""
Author:  Cory Padgett
Advisor: Dr. Jeffrey Fung
Email:   cpadge4@clemson.edu
"""

#imports
import functions as f
import numpy as np
import matplotlib.pyplot as plt
import read_penguin as r
import functions_plots as fp


#dir creation, indexing frames, and setting run info
#------------------------------------------------------------------------------
#edit for run starts here
#------------------------------------------------------------------------------
RUN_NAME      = "h35_2p636E_e30"
ANALYSIS_NAME = str(f.T0)+"_"+str(f.ALPHA)+"_"+str(f.N0)+"_"+str(f.BETA)
RUN_TYPE      = "run_stationary"  #"run_stationary" or "run_rotate"
ORBIT_NUM     = 0                 #Used only in run_rotate run type
FRAME_START   = 0                 #initial frame to pull
FRAME_END     = 500               #final frame -1 to pull
DIM           = '2D'              #"2D" or "3D"
IMAX          = 1200              #x dir grid size
JMAX          = 2160              #y dir grid size
kmax          = 0                 #z dir grid size
LABEL = "h35_2p636E_e30_a-20_b-20_OA_PPM4"   #simulation label
#------------------------------------------------------------------------------
#edit for run ends here
#------------------------------------------------------------------------------


file_name = ['X Position', 'Y Position', 'Density', 'Pressure',
            'X Velocity' ,'Y Velocity', 'Temperatures', 'Multi Plot',
            'Proj-Pol Velocity', 'Polar Density', 'Line Density',
            'Col Dens', 'Optical Depth', 'Pericenter']


f.dir_check(RUN_NAME, ANALYSIS_NAME, file_name, DIM)
vel_j, ein_j, eng_j, g1j, g0j, temp, part = f.load_spec_data()


for k in range(FRAME_START,FRAME_END):
    grids, lengths, data = f.load_PEnGUIn_2d(LABEL, IMAX, JMAX, k)
    radius, theta = grids
    x_len, y_len  = lengths
    time, pressure, peng_density, temp, r_vel, t_vel = data

    #computational constants
    frame = k - FRAME_START
    roll_num = 12 * frame
    #number to roll array's by for frame num

    radius_1 = np.ones((y_len-1,)) * f.INNER
    radius_2 = np.ones((y_len-1,)) * f.OUTER

    #----------------------------------------------------------------------------------------------
    #Line Desnsity Calculation
    #----------------------------------------------------------------------------------------------
    p_temp  = f.T0 * (radius/f.R_IN)**f.ALPHA
    p_dens  = f.N0 * (radius/f.R_IN)**f.BETA
    temp    = p_temp
    density = p_dens

    if RUN_TYPE=="run_rotate":
        temp  = np.roll(temp,  roll_num, axis=0)
        density = np.roll(density, roll_num, axis=0)

    vzp       = f.incvel_2D(r_vel, t_vel, theta, RUN_TYPE, roll_num)
    part_func = f.part_func_1D(temp, part, y_len)
    flux_dens, optical_t = f.flux_cal(density, vel_j, ein_j, eng_j, g1j, g0j, temp, part_func)
    int_j     = f.int_j_cal(flux_dens, radius, f.k_vel, vzp)
    print("Calculations Done!")

    #----------------------------------------------------------------------------------------------
    #Plotting
    #----------------------------------------------------------------------------------------------
    parameters = RUN_TYPE, RUN_NAME, ANALYSIS_NAME, time, k, frame, ORBIT_NUM
    varibles    = f.k_vel, radius, theta, peng_density, int_j, radius_1, radius_2
    fp.plot_line_obs_dens(parameters, varibles)
    fp.plot_line_obs(parameters, f.k_vel, int_j)

    fp.plot_pol_vzp(parameters, radius, theta, radius_1, radius_2, vzp)
    fp.plot_pol_dens(parameters, radius, theta, radius_1, radius_2, peng_density)

    vel_type = "X Velocity"
    fp.plot_vel(parameters, radius, theta, r_vel, ORBIT_NUM, vel_type)
    vel_type = "Y Velocity"
    fp.plot_vel(parameters, radius, theta, t_vel, ORBIT_NUM, vel_type)
