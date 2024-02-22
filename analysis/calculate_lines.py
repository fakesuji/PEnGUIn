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
RUN_NAME      = "h50_2p623E_e40_a-30"
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
vel_j, ein_j, eng_j, g1j, g0j, tempature, part = f.load_spec_data()

int_obs = np.loadtxt("Data/CITau_obs_line.txt")
new_obs = f.remove_spike(int_obs)

data_full = np.loadtxt("Data/CITau_obs_line_full.txt")
std_data_full = np.append(data_full[0:330], data_full[580:-1])
std_full  = np.std(std_data_full)
print("UN-NORM STD of OBS: ",std_full)

var_full = np.sqrt(std_full)
print("UN-NORM Var of OBS:",var_full)

for k in range(FRAME_START,FRAME_END):
    grids, lengths, data = f.load_PEnGUIn_2d(LABEL, IMAX, JMAX, k)
    radius, theta = grids
    x_len, y_len  = lengths
    time, pressure, peng_density, peng_temp, r_vel, t_vel = data
    
    radius, r_vel, t_vel = f.scale_data(radius, r_vel, t_vel)
    
    if x_len == 0:
        print("Orbit " + str(k) + " does not exist")
        continue

    #computational constants for rotation
    frame = 0
    roll_num = 0

    radius_1 = np.ones((y_len-1,)) * f.INNER
    radius_2 = np.ones((y_len-1,)) * f.OUTER
    
    #----------------------------------------------------------------------------------------------
    #Ecc and Angle Calculations
    #----------------------------------------------------------------------------------------------
    ecc, ang = f.ecc_calculation(theta, r_vel, t_vel)
    disk_misalign, ecc_avg, ang_avg, ang = f.ecc_ang_avgs(radius[0,:], ecc, ang, peng_density)
    
    constants_values = [k * f.orbit_interval, frame, disk_misalign, ecc_avg[0], ecc_avg[1], ang_avg[0], ang_avg[1]]  
    ecc_file_path = "Plots"+"/"+RUN_NAME+"/"+ANALYSIS_NAME+"/"+"ecc_ang.csv"
    f.save_ecc_ang_to_csv(constants_values, ecc_file_path)
    
    angle_1 = np.ones((x_len-1,)) * ang_avg[0]
    angle_2 = np.ones((x_len-1,)) * ang_avg[1]

    
    #----------------------------------------------------------------------------------------------
    #Line Desnsity Calculations
    #----------------------------------------------------------------------------------------------
    p_temp  = f.T0 * (radius/f.RADIUS_IN)**f.ALPHA
    p_dens  = f.N0 * (radius/f.RADIUS_IN)**f.BETA
    temp    = p_temp
    density = p_dens

    if RUN_TYPE=="run_rotate":
        temp    = np.roll(temp,    roll_num, axis=0)
        density = np.roll(density, roll_num, axis=0)
        peng_density = np.roll(peng_density, roll_num, axis=0)
        t_vel_roll   = np.roll(t_vel,  roll_num, axis=0)
        r_vel_roll   = np.roll(r_vel,  roll_num, axis=0)

    vzp       = f.incvel_2D(r_vel, t_vel, theta, RUN_TYPE, roll_num)
    part_func = f.part_func_1D(temp, part, x_len)
    flux_dens = f.flux_cal(density, vel_j, ein_j, eng_j, g1j, g0j, temp, part_func)
    int_j     = f.int_j_cal(flux_dens, radius, f.k_vel, vzp)
    int_j     = int_j / f.trap_area(f.k_vel, int_j) #2d area check to see if it is the same
    csv_file_path = "Plots/" + RUN_NAME + "/" + ANALYSIS_NAME + "/" + 'IJ_Data.csv'
    f.save_line_to_csv(int_j, k, frame, csv_file_path)
    
    new_j = f.remove_spike(int_j)
    new_k = f.remove_spike(f.k_vel)
    norm_obs, norm_j, norm_k, norm1, norm2 = f.equal_area_data(vzp, new_j, new_obs, new_k)
    norm_std = std_full / norm1
    norm_var = norm_std**2
    
    
    #----------------------------------------------------------------------------------------------
    #Stats Calculations
    #----------------------------------------------------------------------------------------------
    chi_2 = np.sum(abs(((norm_j - norm_obs)**(2.0) / norm_var))) / (len(norm_k) - 3)
    L_2   = np.sqrt(np.sum((norm_j - norm_obs)**(2.0)))
    L_inf = np.max(norm_j - norm_obs)

    constants_values = [k*10, frame, norm_std, norm_var, chi_2, L_2, L_inf]  
    csv_file_path = "Plots"+"/"+RUN_NAME+"/"+ANALYSIS_NAME+"/"+"stats.csv"
    f.save_stats_to_csv(constants_values, csv_file_path)


    #----------------------------------------------------------------------------------------------
    #Plotting
    #----------------------------------------------------------------------------------------------
    parameters = RUN_TYPE, RUN_NAME, ANALYSIS_NAME, time, k * f.orbit_interval, frame, ORBIT_NUM
    variables   = radius, theta, peng_density, radius_1, radius_2, angle_1, angle_2
    
    fp.plot_line_obs_dens(parameters, variables, norm_obs, norm_j, norm_k)
    fp.plot_line_obs_std_dens(parameters, variables, norm_obs, norm_j, norm_k, norm_std)
    fp.plot_line_obs(parameters, norm_obs, norm_j, norm_k)
    fp.plot_line_obs_std(parameters, norm_obs, norm_j, norm_k, norm_std)

    fp.plot_pol_vzp(parameters, radius, theta, radius_1, radius_2, vzp)
    fp.plot_pol_dens(parameters, variables)

    if RUN_TYPE=="run_rotate":
        vel_type = "X Velocity"
        fp.plot_pol_vel(parameters, radius, theta, r_vel_roll, vel_type)
        vel_type = "Y Velocity"
        fp.plot_pol_vel(parameters, radius, theta, t_vel_roll, vel_type)

    if RUN_TYPE=="run_stationary":
        vel_type = "X Velocity"
        fp.plot_pol_vel(parameters, radius, theta, r_vel, vel_type)
        vel_type = "Y Velocity"
        fp.plot_pol_vel(parameters, radius, theta, t_vel, vel_type)