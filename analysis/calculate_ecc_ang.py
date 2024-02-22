"""
Author:  Cory Padgett
Advisor: Dr. Jeffrey Fung
Email:   cpadge4@clemson.edu
"""

#imports
import functions as f
import numpy as np
import functions_plots as fp
import sys

#dir creation, indexing frames, and setting run info
#------------------------------------------------------------------------------
#edit for run starts here
#------------------------------------------------------------------------------
RUN_NAME      = "h20_2p623E_e80_a-30"
ANALYSIS_NAME = str(f.T0)+"_"+str(f.ALPHA)+"_"+str(f.N0)+"_"+str(f.BETA)
RUN_TYPE      = "run_stationary"  #"run_stationary" or "run_rotate"
ORBIT_NUM     = 0                 #Used only in run_rotate run type
FRAME_START   = 0
FRAME_END     = 301
DIM           = '2D'              #"2D" or "3D"
IMAX          = 1200              #x dir grid size
JMAX          = 2160              #y dir grid size
kmax          = 0                 #z dir grid size
LABEL = RUN_NAME+"_b-20_OA_PPM4"   #simulation label
#------------------------------------------------------------------------------
#edit for run ends here
#------------------------------------------------------------------------------

file_name = ['X Velocity' ,'Y Velocity', 'Multi Plot',
            'Other', 'Polar Density', 'Line Density']

f.dir_check(RUN_NAME, ANALYSIS_NAME, file_name)

for k in range(FRAME_START,FRAME_END):
    grids, lengths, data = f.load_PEnGUIn_2d(LABEL, IMAX, JMAX, k)
    radius, theta = grids
    x_len, y_len  = lengths
    time, pressure, peng_density, peng_temp, r_vel, t_vel = data
    
    radius, r_vel, t_vel = f.scale_data(radius, r_vel, t_vel)
    
    if x_len == 0:
        print("Orbit " + str(k * f.orbit_interval) + " does not exist")
        continue
        
    frame = 0
    radius_1 = np.ones((y_len-1,)) * f.INNER
    radius_2 = np.ones((y_len-1,)) * f.OUTER

    ecc, ang = f.ecc_calculation(theta, r_vel, t_vel)
    disk_misalign, ecc_avg, ang_avg, ang = f.ecc_ang_avgs(radius[0,:], ecc, ang, peng_density)
    
    constants_values = [k * f.orbit_interval, frame, disk_misalign, ecc_avg[0], ecc_avg[1], ang_avg[0], ang_avg[1]]  
    ecc_file_path = "Plots"+"/"+RUN_NAME+"/"+ANALYSIS_NAME+"/"+"ecc_ang.csv"
    f.save_ecc_ang_to_csv(constants_values, ecc_file_path)
    
    angle_1 = np.ones((x_len-1,)) * ang_avg[0]
    angle_2 = np.ones((x_len-1,)) * ang_avg[1]

    #----------------------------------------------------------------------------------------------
    #Plotting
    #----------------------------------------------------------------------------------------------
    parameters = RUN_TYPE, RUN_NAME, ANALYSIS_NAME, time, k * f.orbit_interval, frame, ORBIT_NUM
    variables   = radius, theta, peng_density, radius_1, radius_2, angle_1, angle_2

    fp.plot_dens_ecc_ang(parameters, variables, ecc, ang)
    fp.plot_dens_ecc_ang_full(parameters, variables, ecc, ang)

    fp.plot_pol_dens(parameters, variables)
    fp.plot_pol_dens_full(parameters, variables)

    vel_type = "X Velocity"
    fp.plot_pol_vel(parameters, radius, theta, r_vel, vel_type)
    vel_type = "Y Velocity"
    fp.plot_pol_vel(parameters, radius, theta, t_vel, vel_type)