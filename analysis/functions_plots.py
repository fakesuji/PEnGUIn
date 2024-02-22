"""
Author:  Cory Padgett
Advisor: Dr. Jeffrey Fung
Email:   cpadge4@clemson.edu
"""

#imports
import numpy as np
import functions as f
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


#-------------------------------------------------------------------------------------------------------
#constants
#-------------------------------------------------------------------------------------------------------
#general constants
PLANCK     = 6.626068 * 10**(-27)      #erg s       Planck’s Constant
BOLTZ      = 0.6950348                 #cm^-1 K^-1  Boltzman’s Constant
BOLTZ_CGS  = 1.380649 * 10**(-16)      #erg K^-1  Boltzman’s Constant
C          = 2.997925 * 10**10         #cm s^-1     Speed of light in vacuum
G          = 6.674 * 10**(-8)          #cm^3 g^−1 s^−2
R          = 8.31446261815324 * 10**7  #erg K^-1 mol^-1
MOL        = 6.022 * 10*23
AU_CM      = 1.496 * 10**13     #cm
AU_KM      = 1.496 * 10**8      #km
MASS_SOL   = 1.989 * 10**33     #g
MASS_JUP   = 1.898 * 10**30     #g
H_CO_RATIO = 1000.0
H_MASS     = 2.01568            #g/mol


#system constants
MASS_STAR      = f.MASS_STAR
MASS_PLANET    = f.MASS_PLANET
PLANET_SCALING = f.PLANET_SCALING
INCLIN  = f.INCLIN
T0      = f.T0
ALPHA   = f.ALPHA
N0      = f.N0
BETA    = f.BETA
B       = f.B
B_CGS   = f.B_CGS
INNER   = f.INNER
OUTER   = f.OUTER

#plotting constants
vel_min = -200
vel_max =  200
lin_min = -0.002
lin_max =  0.01
ecc_min =  0.0
ecc_max =  0.25
RADIUS_IN    =  0.25 * f.PLANET_SCALING
RADIUS_OUT   =  10.0 * f.PLANET_SCALING

legend_size    = 12
axis_size      = 12
title_size     = 16
sub_title_size = 14

padding = 2.5
ticks = []
tick_angle = -30.0
reduction_factor = 4 #reduces rmax by 1/#

#arrays
k_vel = np.linspace(-300,300,601)


#-------------------------------------------------------------------------------------------------------
#plotting functions
#-------------------------------------------------------------------------------------------------------
#plot line density spectrum ----------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def plot_line(parameters, int_j):
    #plot ro-vibrational line density
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters

    fig, ax = plt.subplots()
    ax.plot(k_vel,int_j, label="PEnGUIUn")
    plt.legend(fontsize = legend_size)
    plt.xlabel("Velocity [km/s]", fontsize = axis_size)
    plt.ylabel("Normalized Flux", fontsize = axis_size)
    plt.xlim(vel_min, vel_max)
    plt.ylim(lin_min, lin_max)
    plt.grid()
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.title("Line Density at Orbit " + str(k), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.title("Line Density at Orbit " + str(orbit_num) + " and Frame " + str(frame), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_obs_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_obs(parameters, int_obs, int_j, int_k):
    #plot ro-vibrational line density with observation line overlay
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters

    fig, ax = plt.subplots()
    ax.plot(int_k, int_j,   label="PEnGUIUn")
    ax.plot(int_k, int_obs, label="CI Tau Observed")
    plt.legend(fontsize = legend_size)
    plt.xlabel("Velocity [km/s]", fontsize = axis_size)
    plt.ylabel("Normalized Flux", fontsize = axis_size)
    plt.xlim(vel_min, vel_max)
    plt.ylim(lin_min, lin_max)
    plt.grid()
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.title("Line Density at Orbit " + str(k), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_obs_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.title("Line Density at Orbit " + str(orbit_num) + " and Frame " + str(frame), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_obs_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_obs_std(parameters, int_obs, int_j, int_k, norm_std):
    #plot ro-vibrational line density with observation line overlay
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters

    std1_upper = int_obs + norm_std
    std1_lower = int_obs - norm_std
    std2_upper = int_obs + 2*norm_std
    std2_lower = int_obs - 2*norm_std

    fig, ax = plt.subplots()
    ax.plot(int_k, int_obs, label="CI Tau Observed", color='orange')
    ax.fill_between(int_k, std1_lower, std1_upper, color='g', alpha=0.6, label="1 STD")
    ax.fill_between(int_k, std2_lower, std2_upper, color='g', alpha=0.3, label="2 STD")
    ax.plot(int_k, int_j, color='b', label="PEnGUIUn")
    
    plt.legend(fontsize = legend_size)
    plt.xlim(vel_min, vel_max)
    plt.ylim(lin_min, lin_max)
    plt.grid()
    plt.xlabel("Velocity [km/s]", fontsize = axis_size)
    plt.ylabel("Normalized Flux", fontsize = axis_size)
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.title("Line Density at Orbit " + str(k), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_obs_std_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.title("Line Density at Orbit " + str(orbit_num) + " and Frame " + str(frame), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_obs_std_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_seg(parameters, int_j, int_k, int_b):
    #plot ro-vibrational line density and in spatial segments
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    
    b_len = len(int_b[:])
    b_label = []

    for i in range(b_len):
        b_label.append("Segment " + str(i))

    fig, ax =plt.subplots()
    ax.plot(int_k, int_j,        label="Sum of Segments")
    for i in range(b_len):
        ax.plot(int_k, int_b[i], label=b_label[i])
    
    plt.legend(fontsize = legend_size)
    plt.xlim(vel_min, vel_max)
    plt.ylim(lin_min, lin_max)
    plt.grid()
    plt.xlabel("Velocity [km/s]", fontsize = axis_size)
    plt.ylabel("Normalized Flux", fontsize = axis_size)
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.title("Line Density at Orbit " + str(k), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_seg_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.title("Line Density at Orbit " + str(orbit_num) + " and Frame " + str(frame), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_seg_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_obs_seg(parameters, int_obs, int_j, int_k, int_b):
    #plot ro-vibrational line density with observation line overlay and in spatial segments
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    
    b_len = len(int_b[:])
    b_label = []

    for i in range(b_len):
        b_label.append("Segment " + str(i))

    fig, ax =plt.subplots()
    ax.plot(int_k, int_j,   label="Sum of Segments")
    ax.plot(int_k, int_obs, label="CI Tau Observed")
    for i in range(b_len):
        ax.plot(int_k, int_b[i], label=b_label[i])
    
    plt.legend(fontsize = legend_size)
    plt.xlim(vel_min, vel_max)
    plt.ylim(lin_min, lin_max)
    plt.grid()
    plt.xlabel("Velocity [km/s]", fontsize = axis_size)
    plt.ylabel("Normalized Flux", fontsize = axis_size)
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.title("Line Density at Orbit " + str(k), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_obs_seg_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.title("Line Density at Orbit " + str(orbit_num) + " and Frame " + str(frame), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_obs_seg_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


#-------------------------------------------------------------------------------------------------------
#polar plots -------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def plot_pol_vel(parameters, radius, theta, vel, vel_type):
    #Polar velocity plot
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    
    if vel_type == "X Velocity":
        v_type = "x"
        vel_min = -15
        vel_max =  15
        
    if vel_type == "Y Velocity":
        v_type = "y"
        vel  = vel - np.sqrt(G * MASS_STAR / (radius * AU_CM))/100000
        vel_min = -15
        vel_max =  15

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6,5))
    ax.grid(False)
    plt.pcolormesh(theta, radius, vel, vmin = vel_min, vmax = vel_max)
    
    ax.set_rmax(RADIUS_OUT/reduction_factor)
    ax.set_rticks(ticks)
    ax.set_rlabel_position(tick_angle)
    plt.colorbar(label="[km/s]")
    plt.tight_layout(pad=padding)
    
    
    if run_type == "run_stationary":
        ax.set_title("PEnGUIn " + vel_type + " Velocities at Orbit " + str(k), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + vel_type + "/" + v_type + "_vel_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax.set_title("PEnGUIn " + vel_type + " Velocities at Orbit " + str(orbit_num) + " and Frame " + str(frame), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + vel_type + "/" + v_type + "_vel_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_pol_dens(parameters, variables):
    #Polar Density Map
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6,5))
    ax.grid(False)
    plt.pcolormesh(theta, radius, np.log10(density), vmin=-2, vmax=2)
    
    ax.plot(angle_1, radius[0,:], color='k', linestyle='dashed', label="Inner Peri")
    ax.plot(angle_2, radius[0,:], color='r', linestyle='dashed', label="Outer Peri")
    
    ax.set_rmax(RADIUS_OUT/reduction_factor)
    ax.set_rticks(ticks)
    ax.set_rlabel_position(tick_angle)
    
    cbar = plt.colorbar(ax=ax, label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    plt.legend(bbox_to_anchor=(-0.25, 1.05), loc='upper left', fontsize = legend_size)
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        ax.set_title("Calculated PEnGUIn Densities at Orbit " + str(k), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Polar Density" + "/pol_dens_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax.set_title("Calculated PEnGUIn Densities at Orbit " + str(orbit_num) + " and Frame " + str(frame), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Polar Density" + "/pol_dens_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return

def plot_pol_dens_full(parameters, variables):
    #Polar Density Map
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6,5))
    ax.grid(False)
    plt.pcolormesh(theta, radius, np.log10(density), vmin=-2, vmax=2)
    
    ax.plot(angle_1, radius[0,:], color='k', linestyle='dashed', label="Inner Peri")
    ax.plot(angle_2, radius[0,:], color='r', linestyle='dashed', label="Outer Peri")
    
    ax.set_rmax(RADIUS_OUT)
    ax.set_rticks(ticks)
    ax.set_rlabel_position(tick_angle)
    
    cbar = plt.colorbar(ax=ax, label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    plt.legend(bbox_to_anchor=(-0.25, 1.05), loc='upper left', fontsize = legend_size)
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        ax.set_title("Calculated PEnGUIn Densities at Orbit " + str(k), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Polar Density" + "/pol_dens_full_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax.set_title("Calculated PEnGUIn Densities at Orbit " + str(orbit_num) + " and Frame " + str(frame), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Polar Density" + "/pol_dens_full_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_pol_vzp(parameters, radius, theta, radius_1, radius_2, vzp):
    #Polar Inclined Velocity Map
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.grid(False)
    plt.pcolormesh(theta, radius, vzp, vmin=-150.0, vmax=150.0)
    
    ax.set_rmax(RADIUS_OUT/reduction_factor)
    ax.set_rticks(ticks)
    ax.set_rlabel_position(tick_angle)
    
    plt.colorbar(label="[km/s]")
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        ax.set_title("Projected PEnGUIn Velocities at Orbit " + str(k), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Other" + "/pol_vzp_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax.set_title("Projected PEnGUIn Velocities at Orbit " + str(orbit_num) + " and Frame " + str(frame), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Other" + "/pol_vzp_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_pol_temp(parameters, radius, theta, radius_1, radius_2, temp):
    #Polar Temperature Map
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.grid(False)
    plt.pcolormesh(theta, radius, np.log10(temp))
    
    ax.set_rmax(RADIUS_OUT/reduction_factor)
    ax.set_rticks(ticks)
    ax.set_rlabel_position(tick_angle)
    
    plt.colorbar(label="[$\log_{10}$(k)]")
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        ax.set_title("Surface Temperatures at Orbit " + str(k), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Other" + "/temp_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax.set_title("Surface Temperatures at Orbit " + str(orbit_num) + " and Frame " + str(frame), fontsize = title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Other" + "/temp_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


#-------------------------------------------------------------------------------------------------------
#Multi-info Plots --------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def plot_line_obs_dens_ecc(parameters, variables, ecc, int_obs, int_j, int_k):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables

    fig, _axs = plt.subplots(nrows=1, ncols=3, figsize=(18, 5))
    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)
    axs = _axs.flatten() 

    cset0 = axs[0].plot(int_k, int_j,   label="PEnGUIn")
    axs[0].plot(int_k,         int_obs, label="CI Tau Observed")
    cset1 = axs[1].pcolormesh(radius, theta, np.log10(density), vmin=-2, vmax=2)
    cset2 = axs[2].scatter(radius[0,:], ecc, s=1, c='r')

    axs[1].plot(radius[0,:], angle_1, color='k', linestyle='dashed', label="Inner Peri")
    axs[1].plot(radius[0,:], angle_2, color='r', linestyle='dashed', label="Outer Peri")
    
    axs[0].set_xlabel("Velocities [km/s]", fontsize = axis_size)
    axs[0].set_ylabel("Normalized Flux", fontsize = axis_size)
    axs[0].legend(fontsize = legend_size)
    axs[0].set_xlim(vel_min, vel_max)
    axs[0].set_ylim(lin_min, lin_max)
    axs[0].grid()

    axs[1].set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)
    axs[1].set_xlabel("Radius [au]", fontsize = axis_size)
    axs[1].set_xlim(RADIUS_IN, RADIUS_OUT)
    axs[1].set_ylabel("Theta [rads]", fontsize = axis_size)

    axs[2].set_title("PEnGUIn Disk Ecc", fontsize = sub_title_size)
    axs[2].set_xlabel("Radius [au]", fontsize = axis_size)
    axs[2].set_xlim(RADIUS_IN, RADIUS_OUT)
    axs[2].set_ylim(ecc_min, ecc_max)
    axs[2].set_ylabel("Ecc", fontsize = axis_size)

    #axs[1].plot(radius_1, theta[:,0], color='k', linestyle='dotted')
    #axs[1].plot(radius_2, theta[:,0], color='k', linestyle='dotted')

    #axs[2].vlines(INNER, ecc_min, ecc_max, color='k', linestyle='dotted')
    #axs[2].vlines(OUTER, ecc_min, ecc_max, color='k', linestyle='dotted')

    cbar = plt.colorbar(cset1, ax=axs[1], label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    axs[1].legend(bbox_to_anchor=(0.5, 0.95), loc='upper left', fontsize = legend_size)
    
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_ecc_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_ecc_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_obs_std_dens_ecc(parameters, variables, ecc, int_obs, int_j, int_k, norm_std):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables
    
    std1_upper = int_obs + norm_std
    std1_lower = int_obs - norm_std
    std2_upper = int_obs + 2*norm_std
    std2_lower = int_obs - 2*norm_std

    fig, _axs = plt.subplots(nrows=1, ncols=3, figsize=(18, 5))
    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)
    axs = _axs.flatten() 

    cset0 = axs[0].plot(int_k, int_j,   label="PEnGUIn")
    axs[0].plot(int_k,         int_obs, label="CI Tau Observed")
    axs[0].fill_between(int_k, std1_lower, std1_upper, color='g', alpha=0.6, label="1 STD")
    axs[0].fill_between(int_k, std2_lower, std2_upper, color='g', alpha=0.3, label="2 STD")
    cset1 = axs[1].pcolormesh(radius, theta, np.log10(density), vmin=-2, vmax=2)
    cset2 = axs[2].scatter(radius[0,:], ecc, s=1, c='r')

    axs[0].set_xlabel("Velocities [km/s]", fontsize = axis_size)
    axs[0].set_ylabel("Normalized Flux", fontsize = axis_size)
    axs[0].set_xlim(vel_min, vel_max)
    axs[0].set_ylim(lin_min, lin_max)
    axs[0].grid()

    axs[1].plot(radius[0,:], angle_1, color='k', linestyle='dashed', label="Inner Peri")
    axs[1].plot(radius[0,:], angle_2, color='r', linestyle='dashed', label="Outer Peri")
    
    axs[1].set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)
    axs[1].set_xlabel("Radius [au]", fontsize = axis_size)
    axs[1].set_xlim(RADIUS_IN, RADIUS_OUT)
    axs[1].set_ylabel("Theta [rads]", fontsize = axis_size)

    axs[2].set_title("PEnGUIn Disk Ecc", fontsize = sub_title_size)
    axs[2].set_xlabel("Radius [au]", fontsize = axis_size)
    axs[2].set_xlim(RADIUS_IN, RADIUS_OUT)
    axs[2].set_ylim(ecc_min, ecc_max)
    axs[2].set_ylabel("Ecc", fontsize = axis_size)

    axs[1].plot(radius_1, theta[:,0], color='k', linestyle='dotted')
    axs[1].plot(radius_2, theta[:,0], color='k', linestyle='dotted')

    #axs[2].vlines(INNER, ecc_min, ecc_max, color='k', linestyle='dotted')
    #axs[2].vlines(OUTER, ecc_min, ecc_max, color='k', linestyle='dotted')
    
    axs[0].legend(fontsize = legend_size)

    cbar = plt.colorbar(cset1, ax=axs[1], label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    axs[1].legend(bbox_to_anchor=(0.5, 0.95), loc='upper left', fontsize = legend_size)
    
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_std_ecc_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_std_ecc_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return
    

def plot_line_obs_dens_ecc_ang(parameters, variables, ecc, ang, int_obs, int_j, int_k):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables
    
    fig = plt.figure(figsize=(10, 9))

    gs = GridSpec(7, 6, figure=fig)
    ax1 = fig.add_subplot(gs[0:3, 0:3])
    ax2 = fig.add_subplot(gs[0:3, 3:6])
    ax3 = fig.add_subplot(gs[3:5, :])
    ax4 = fig.add_subplot(gs[5:7, :])


    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)

    ax1.plot(int_k, int_j,   label="PEnGUIn")
    ax1.plot(int_k, int_obs, label="CI Tau Observed")

    cset2 = ax2.pcolormesh(radius, theta, np.log10(density), vmin=-2, vmax=2)
    ax2.plot(radius_1, theta[:,0], color='k', linestyle='dotted')     #plot visual inner gap boundary
    ax2.plot(radius_2, theta[:,0], color='k', linestyle='dotted')     #plot visual outer gap boundary

    ax1.set_xlabel("Velocities [km/s]", fontsize = axis_size)
    ax1.set_ylabel("Normalized Flux", fontsize = axis_size)
    ax1.set_xlim(vel_min, vel_max)
    ax1.set_ylim(lin_min, lin_max)
    ax1.grid()

    ax2.plot(radius[0,:], angle_1, color='k', linestyle='dashed', label="Inner Peri")
    ax2.plot(radius[0,:], angle_2, color='r', linestyle='dashed', label="Outer Peri")
    
    ax2.set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)
    ax2.set_xlabel("Radius [au]", fontsize = axis_size)
    ax2.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax2.set_ylabel("Theta [rads]", fontsize = axis_size)

    ax3.scatter(radius[0, :], ecc, s=1, c='r')
    ax3.set_xlabel("Radius [au]", fontsize = axis_size)
    ax3.set_ylabel("Ecc", fontsize = axis_size)
    ax3.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax3.set_ylim(ecc_min, ecc_max)
    ax3.set_title("PEnGUIn Disk Ecc and and Arguement of Periapsis", fontsize = sub_title_size)
    #ax3.vlines(INNER, ecc_min, ecc_max, color='k', linestyle='dotted')
    #ax3.vlines(OUTER, ecc_min, ecc_max, color='k', linestyle='dotted')

    ax4.scatter(radius[0, :], ang, s=1, c='k')
    ax4.set_xlabel("Radius [au]", fontsize = axis_size)
    ax4.set_ylabel("$\\omega$ [Rads]", fontsize = axis_size)
    ax4.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax4.set_ylim(0, 2*np.pi)
    #ax4.vlines(INNER, -np.pi, np.pi, color='k', linestyle='dotted')
    #ax4.vlines(OUTER, -np.pi, np.pi, color='k', linestyle='dotted')
    
    ax1.legend(fontsize = legend_size)

    cbar = plt.colorbar(cset2, ax=ax2, label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    ax2.legend(bbox_to_anchor=(-0.3, 1.05), loc='upper left', fontsize = legend_size)
    
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_ecc_ang_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_ecc_ang_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return
    
    
def plot_line_obs_std_dens_ecc_ang(parameters, variables, ecc, ang, int_obs, int_j, int_k, norm_std):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables
    
    std1_upper = int_obs + norm_std
    std1_lower = int_obs - norm_std
    std2_upper = int_obs + 2*norm_std
    std2_lower = int_obs - 2*norm_std
    
    fig = plt.figure(figsize=(10, 9))

    gs = GridSpec(7, 6, figure=fig)
    ax1 = fig.add_subplot(gs[0:3, 0:3])
    ax2 = fig.add_subplot(gs[0:3, 3:6])
    ax3 = fig.add_subplot(gs[3:5, :])
    ax4 = fig.add_subplot(gs[5:7, :])


    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)

    ax1.plot(int_k, int_j,   label="PEnGUIn")
    ax1.plot(int_k, int_obs, label="CI Tau Observed")
    ax1.fill_between(int_k, std1_lower, std1_upper, color='g', alpha=0.6, label="1 STD")
    ax1.fill_between(int_k, std2_lower, std2_upper, color='g', alpha=0.3, label="2 STD")

    cset2 = ax2.pcolormesh(radius, theta, np.log10(density), vmin=-2, vmax=2)
    ax2.plot(radius_1, theta[:,0], color='k', linestyle='dotted')     #plot visual inner gap boundary
    ax2.plot(radius_2, theta[:,0], color='k', linestyle='dotted')     #plot visual outer gap boundary

    ax1.set_xlabel("Velocities [km/s]", fontsize = axis_size)
    ax1.set_ylabel("Normalized Flux", fontsize = axis_size)
    ax1.legend(fontsize = legend_size)
    ax1.set_xlim(vel_min, vel_max)
    ax1.set_ylim(lin_min, lin_max)
    ax1.grid()

    ax2.plot(angle_1, radius[0,:], color='k', linestyle='dashed', label="Inner Peri")
    ax2.plot(angle_2, radius[0,:], color='r', linestyle='dashed', label="Outer Peri")
    
    ax2.set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)
    ax2.set_xlabel("Radius [au]", fontsize = axis_size)
    ax2.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax2.set_ylabel("Theta [rads]", fontsize = axis_size)

    ax3.scatter(radius[0, :], ecc, s=1, c='r')
    ax3.set_xlabel("Radius [au]", fontsize = axis_size)
    ax3.set_ylabel("Ecc", fontsize = axis_size)
    ax3.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax3.set_ylim(ecc_min, ecc_max)
    ax3.set_title("PEnGUIn Disk Ecc and and Arguement of Periapsis", fontsize = sub_title_size)
    #ax3.vlines(INNER, ecc_min, ecc_max, color='k', linestyle='dotted')
    #ax3.vlines(OUTER, ecc_min, ecc_max, color='k', linestyle='dotted')

    ax4.scatter(radius[0, :], ang, s=1, c='k')
    ax4.set_xlabel("Radius [au]", fontsize = axis_size)
    ax4.set_ylabel("$\\omega$ [Rads]", fontsize = axis_size)
    ax4.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax4.set_ylim(0, 2*np.pi)
    #ax4.vlines(INNER, -np.pi, np.pi, color='k', linestyle='dotted')
    #ax4.vlines(OUTER, -np.pi, np.pi, color='k', linestyle='dotted')
    
    ax1.legend(fontsize = legend_size)

    cbar = plt.colorbar(cset2, ax=ax2, label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    ax2.legend(bbox_to_anchor=(0.5, 0.95), loc='upper left', fontsize = legend_size)
    
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_std_dens_ecc_ang_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_std_dens_ecc_ang_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_dens_ecc_ang(parameters, variables, ecc, ang):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables
    
    fig = plt.figure(figsize=(16, 5))
    gs = GridSpec(2, 4, figure=fig)
    ax1 = fig.add_subplot(gs[0:2, 0:2], projection='polar')
    ax2 = fig.add_subplot(gs[0, 2:4])
    ax3 = fig.add_subplot(gs[1, 2:4])

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)
    ax1.grid(False)

    cset1 = ax1.pcolormesh(theta, radius, np.log10(density), vmin=-2, vmax=2)
    ax1.plot(angle_1, radius[0,:], color='k', linestyle='dashed', label="Inner Peri")
    ax1.plot(angle_2, radius[0,:], color='r', linestyle='dashed', label="Outer Peri")
    
    ax1.set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)
    ax1.set_rticks([0.1, 0.2, 0.3])
    ax1.set_rmax(RADIUS_OUT/reduction_factor)

    ax2.scatter(radius[0, :], ecc, s=1, c='r')
    ax2.set_xlabel("Radius [au]", fontsize = axis_size)
    ax2.set_ylabel("Ecc", fontsize = axis_size)
    ax2.set_xlim(RADIUS_IN, RADIUS_OUT/reduction_factor)
    ax2.set_ylim(ecc_min, ecc_max)
    ax2.set_title("PEnGUIn Disk Ecc and and Arguement of Periapsis", fontsize = sub_title_size)
    #ax2.vlines(INNER, ecc_min, ecc_max, color='k', linestyle='dotted')
    #ax2.vlines(OUTER, ecc_min, ecc_max, color='k', linestyle='dotted')

    ax3.scatter(radius[0, :], ang, s=1, c='k')
    ax3.set_xlabel("Radius [au]", fontsize = axis_size)
    ax3.set_ylabel("$\\omega$ [Rads]", fontsize = axis_size)
    ax3.set_xlim(RADIUS_IN, RADIUS_OUT/reduction_factor)
    ax3.set_ylim(0, 2*np.pi)
    #ax3.vlines(INNER, -np.pi, np.pi, color='k', linestyle='dotted')
    #ax3.vlines(OUTER, -np.pi, np.pi, color='k', linestyle='dotted')

    cbar = plt.colorbar(cset1, ax=ax1, label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    ax1.legend(bbox_to_anchor=(-0.3, 1.05), loc='upper left', fontsize = legend_size)
    
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/dens_ecc_ang_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/dens_ecc_ang_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_dens_ecc_ang_full(parameters, variables, ecc, ang):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables
    
    fig = plt.figure(figsize=(16, 5))
    gs = GridSpec(2, 4, figure=fig)
    ax1 = fig.add_subplot(gs[0:2, 0:2], projection='polar')
    ax2 = fig.add_subplot(gs[0, 2:4])
    ax3 = fig.add_subplot(gs[1, 2:4])

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)
    ax1.grid(False)

    cset1 = ax1.pcolormesh(theta, radius, np.log10(density), vmin=-2, vmax=2)
    ax1.plot(angle_1, radius[0,:], color='k', linestyle='dashed', label="Inner Peri")
    ax1.plot(angle_2, radius[0,:], color='r', linestyle='dashed', label="Outer Peri")
    #ax1.plot(theta[:,0], radius_1, color='k', linestyle='dotted')
    #ax1.plot(theta[:,0], radius_2, color='k', linestyle='dotted')
    
    ax1.set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)
    ax1.set_rticks(ticks)
    ax1.set_rmax(RADIUS_OUT)

    ax2.scatter(radius[0, :], ecc, s=1, c='r')
    ax2.set_xlabel("Radius [au]", fontsize = axis_size)
    ax2.set_ylabel("Ecc", fontsize = axis_size)
    ax2.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax2.set_ylim(ecc_min, ecc_max)
    ax2.set_title("PEnGUIn Disk Ecc and and Arguement of Periapsis", fontsize = sub_title_size)
    #ax2.vlines(INNER, ecc_min, ecc_max, color='k', linestyle='dotted')
    #ax2.vlines(OUTER, ecc_min, ecc_max, color='k', linestyle='dotted')

    ax3.scatter(radius[0, :], ang, s=1, c='k')
    ax3.set_xlabel("Radius [au]", fontsize = axis_size)
    ax3.set_ylabel("$\\omega$ [Rads]", fontsize = axis_size)
    ax3.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax3.set_ylim(0, 2*np.pi)
    #ax3.vlines(INNER, -np.pi, np.pi, color='k', linestyle='dotted')
    #ax3.vlines(OUTER, -np.pi, np.pi, color='k', linestyle='dotted')

    cbar = plt.colorbar(cset1, ax=ax1, label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    ax1.legend(bbox_to_anchor=(-0.3, 1.05), loc='upper left', fontsize = legend_size)
    
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/dens_ecc_ang_full_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/dens_ecc_ang_full_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_obs_ecc_ang(parameters, variables, ecc, ang, int_obs, int_j, int_k):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables
    
    fig = plt.figure(figsize=(12, 5))
    gs = GridSpec(2, 4, figure=fig)
    ax1 = fig.add_subplot(gs[0:2, 0:2])
    ax2 = fig.add_subplot(gs[0, 2:4])
    ax3 = fig.add_subplot(gs[1, 2:4])

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)

    ax1.plot(int_k, int_j,   label="PEnGUIn")
    ax1.plot(int_k, int_obs, label="CI Tau Observed")
    ax1.set_xlabel("Velocities [km/s]", fontsize = axis_size)
    ax1.set_ylabel("Normalized Flux", fontsize = axis_size)
    ax1.legend(fontsize = legend_size)
    ax1.set_xlim(vel_min, vel_max)
    ax1.set_ylim(lin_min, lin_max)
    ax1.grid()

    ax2.scatter(radius[0, :], ecc, s=1, c='r')
    ax2.set_xlabel("Radius [au]", fontsize = axis_size)
    ax2.set_ylabel("Ecc", fontsize = axis_size)
    ax2.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax2.set_ylim(ecc_min, ecc_max)
    ax2.set_title("PEnGUIn Disk Ecc and and Arguement of Periapsis", fontsize = sub_title_size)
    #ax2.vlines(INNER, ecc_min, ecc_max, color='k', linestyle='dotted')
    #ax2.vlines(OUTER, ecc_min, ecc_max, color='k', linestyle='dotted')

    ax3.scatter(radius[0, :], ang, s=1, c='k')
    ax3.set_xlabel("Radius [au]", fontsize = axis_size)
    ax3.set_ylabel("$\\omega$ [Rads]", fontsize = axis_size)
    ax3.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax3.set_ylim(0, 2*np.pi)
    #ax3.vlines(INNER, -np.pi, np.pi, color='k', linestyle='dotted')
    #ax3.vlines(OUTER, -np.pi, np.pi, color='k', linestyle='dotted')

    plt.tight_layout(pad=padding)

    if run_type == "run_stationary":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_ecc_ang_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_ecc_ang_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return
    
    
def plot_line_obs_std_ecc_ang(parameters, variables, ecc, ang, int_obs, int_j, int_k, norm_std):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables
    
    std1_upper = int_obs + norm_std
    std1_lower = int_obs - norm_std
    std2_upper = int_obs + 2*norm_std
    std2_lower = int_obs - 2*norm_std
    
    fig = plt.figure(figsize=(12, 5))
    gs = GridSpec(2, 4, figure=fig)
    ax1 = fig.add_subplot(gs[0:2, 0:2])
    ax2 = fig.add_subplot(gs[0, 2:4])
    ax3 = fig.add_subplot(gs[1, 2:4])


    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)

    ax1.plot(int_k, int_j,   label="PEnGUIn")
    ax1.plot(int_k, int_obs, label="CI Tau Observed")
    ax1.fill_between(int_k, std1_lower, std1_upper, color='g', alpha=0.6, label="1 STD")
    ax1.fill_between(int_k, std2_lower, std2_upper, color='g', alpha=0.3, label="2 STD")
    ax1.set_xlabel("Velocities [km/s]", fontsize = axis_size)
    ax1.set_ylabel("Normalized Flux", fontsize = axis_size)
    ax1.legend(fontsize = legend_size)
    ax1.set_xlim(vel_min, vel_max)
    ax1.set_ylim(lin_min, lin_max)
    ax1.grid()

    ax2.scatter(radius[0, :], ecc, s=1, c='r')
    ax2.set_xlabel("Radius [au]", fontsize = axis_size)
    ax2.set_ylabel("Ecc", fontsize = axis_size)
    ax2.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax2.set_ylim(ecc_min, ecc_max)
    ax2.set_title("PEnGUIn Disk Ecc and and Arguement of Periapsis", fontsize = sub_title_size)
    #ax2.vlines(INNER, ecc_min, ecc_max, color='k', linestyle='dotted')
    #ax2.vlines(OUTER, ecc_min, ecc_max, color='k', linestyle='dotted')

    ax3.scatter(radius[0, :], ang, s=1, c='k')
    ax3.set_xlabel("Radius [au]", fontsize = axis_size)
    ax3.set_ylabel("$\\omega$ [Rads]", fontsize = axis_size)
    ax3.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax3.set_ylim(0, 2*np.pi)
    #ax3.vlines(INNER, -np.pi, np.pi, color='k', linestyle='dotted')
    #ax3.vlines(OUTER, -np.pi, np.pi, color='k', linestyle='dotted')

    plt.tight_layout(pad=padding)

    if run_type == "run_stationary":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_std_ecc_ang_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_std_ecc_ang_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_obs_dens_ecc_j_ang(parameters, variables, ecc, ang, int_obs, int_j, int_k):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables
    
    fig = plt.figure(figsize=(18, 5))

    gs = GridSpec(1, 3, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = ax3.twinx()

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)

    ax1.plot(int_k, int_j,   label="PEnGUIn")
    ax1.plot(int_k, int_obs, label="CI Tau Observed")

    cset2 = ax2.pcolormesh(radius, theta, np.log10(density), vmin=-2, vmax=2)
    ax2.plot(radius[0,:], angle_1, color='k', linestyle='dashed', label="Inner Peri")
    ax2.plot(radius[0,:], angle_2, color='r', linestyle='dashed', label="Outer Peri")
    #ax2.plot(radius_1, theta[:,0], color='k', linestyle='dotted')
    #ax2.plot(radius_2, theta[:,0], color='k', linestyle='dotted')

    ax1.set_xlabel("Velocities [km/s]", fontsize = axis_size)
    ax1.set_ylabel("Normalized Flux", fontsize = axis_size)
    ax1.legend(fontsize = legend_size)
    ax1.set_xlim(vel_min, vel_max)
    ax1.set_ylim(lin_min, lin_max)
    ax1.grid()

    ax2.set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)
    ax2.set_xlabel("Radius [au]", fontsize = axis_size)
    ax2.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax2.set_ylabel("Theta [rads]", fontsize = axis_size)

    ax3.scatter(radius[0, :], ecc, s=1, c='r', label="Ecc")
    ax3.set_xlabel("Radius [au]", fontsize = axis_size)
    ax3.set_ylabel("Ecc", fontsize = axis_size)
    ax3.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax3.set_ylim(ecc_min, ecc_max)
    ax3.set_title("PEnGUIn Disk Ecc and Arguement of Periapsis", fontsize = sub_title_size)
    #ax3.vlines(INNER, ecc_min, ecc_max, color='k', linestyle='dotted')
    #ax3.vlines(OUTER, ecc_min, ecc_max, color='k', linestyle='dotted')

    ax4.scatter(radius[0, :], ang, s=1, c='k', label="$\\omega$")
    ax4.set_ylabel("$\\omega$ [Rads]", fontsize = axis_size)
    ax4.set_ylim(0, 2*np.pi)
    
    ax1.legend(fontsize = legend_size)

    lines3, labels3 = ax3.get_legend_handles_labels()
    lines4, labels4 = ax4.get_legend_handles_labels()
    ax3.legend(lines3 + lines4, labels3 + labels4, loc='upper right', fontsize = legend_size)

    cbar = plt.colorbar(cset2, ax=ax2, label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    ax2.legend(bbox_to_anchor=(-0.3, 1.05), loc='upper left', fontsize = legend_size)
    
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_ecc_j_ang_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_ecc_j_ang_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_obs_std_dens_ecc_j_ang(parameters, variables, ecc, ang, int_obs, int_j, int_k, norm_std):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables
    
    std1_upper = int_obs + norm_std
    std1_lower = int_obs - norm_std
    std2_upper = int_obs + 2*norm_std
    std2_lower = int_obs - 2*norm_std
    
    fig = plt.figure(figsize=(18, 5))

    gs = GridSpec(1, 3, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = ax3.twinx()

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)

    ax1.plot(int_k, int_j,   label="PEnGUIn")
    ax1.plot(int_k, int_obs, label="CI Tau Observed")
    ax1.fill_between(int_k, std1_lower, std1_upper, color='g', alpha=0.6, label="1 STD")
    ax1.fill_between(int_k, std2_lower, std2_upper, color='g', alpha=0.3, label="2 STD")

    cset2 = ax2.pcolormesh(radius, theta, np.log10(density), vmin=-2, vmax=2)
    ax2.plot(radius[0,:], angle_1, color='k', linestyle='dashed', label="Inner Peri")
    ax2.plot(radius[0,:], angle_2, color='r', linestyle='dashed', label="Outer Peri")
    #ax2.plot(radius_1, theta[:,0], color='k', linestyle='dotted')
    #ax2.plot(radius_2, theta[:,0], color='k', linestyle='dotted')

    ax1.set_xlabel("Velocities [km/s]", fontsize = axis_size)
    ax1.set_ylabel("Normalized Flux", fontsize = axis_size)
    ax1.legend(fontsize = legend_size)
    ax1.set_xlim(vel_min, vel_max)
    ax1.set_ylim(lin_min, lin_max)
    ax1.grid()

    ax2.set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)
    ax2.set_xlabel("Radius [au]", fontsize = axis_size)
    ax2.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax2.set_ylabel("Theta [rads]", fontsize = axis_size)

    ax3.scatter(radius[0, :], ecc, s=1, c='r', label="Ecc")
    ax3.set_xlabel("Radius [au]", fontsize = axis_size)
    ax3.set_ylabel("Ecc", fontsize = axis_size)
    ax3.set_xlim(RADIUS_IN, RADIUS_OUT)
    ax3.set_ylim(ecc_min, ecc_max)
    ax3.set_title("PEnGUIn Disk Ecc and Arguement of Periapsis", fontsize = sub_title_size)
    #ax3.vlines(INNER, ecc_min, ecc_max, color='k', linestyle='dotted')
    #ax3.vlines(OUTER, ecc_min, ecc_max, color='k', linestyle='dotted')

    ax4.scatter(radius[0, :], ang, s=1, c='k', label="$\\omega$")
    ax4.set_ylabel("$\\omega$ [Rads]", fontsize = axis_size)
    ax4.set_ylim(0, 2*np.pi)
    
    ax1.legend(fontsize = legend_size)

    lines3, labels3 = ax3.get_legend_handles_labels()
    lines4, labels4 = ax4.get_legend_handles_labels()
    ax3.legend(lines3 + lines4, labels3 + labels4, loc='upper right', fontsize = legend_size)

    cbar = plt.colorbar(cset2, ax=ax2, label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    ax2.legend(bbox_to_anchor=(0.5, 0.95), loc='upper left', fontsize = legend_size)
    
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_std_dens_ecc_j_ang_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_std_dens_ecc_j_ang_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_dens_ecc_j_ang(parameters, variables, ecc, ang):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables
    
    fig = plt.figure(figsize=(10.5, 5))

    ax1 = plt.subplot(121, projection='polar')
    ax2 = plt.subplot(122)
    ax3 = ax2.twinx()

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)

    ax1.grid(False)

    cset1 = ax1.pcolormesh(theta, radius, np.log10(density), vmin=-2, vmax=2)
    ax1.plot(angle_1, radius[0,:], color='k', linestyle='dashed', label="Inner Peri")
    ax1.plot(angle_2, radius[0,:], color='r', linestyle='dashed', label="Outer Peri")
    #ax1.plot(theta[:,0], radius_1, color='k', linestyle='dotted')
    #ax1.plot(theta[:,0], radius_2, color='k', linestyle='dotted')
    ax1.set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)
    ax1.set_rmax(RADIUS_OUT/reduction_factor)

    ax2.scatter(radius[0, :], ecc, s=1, c='r')
    ax2.set_xlabel("Radius [au]", fontsize = axis_size)
    ax2.set_ylabel("Ecc", fontsize = axis_size)
    ax2.set_xlim(RADIUS_IN, RADIUS_OUT/reduction_factor)
    ax2.set_ylim(ecc_min, ecc_max)
    ax2.set_title("PEnGUIn Disk Ecc and and Arguement of Periapsis", fontsize = sub_title_size)
    #ax2.vlines(INNER, ecc_min, ecc_max, color='k', linestyle='dotted')
    #ax2.vlines(OUTER, ecc_min, ecc_max, color='k', linestyle='dotted')

    ax3.scatter(radius[0, :], ang, s=1, c='k')
    ax3.set_ylabel("$\\omega$ [Rads]", fontsize = axis_size)
    ax3.set_ylim(0, 2*np.pi)

    cbar = plt.colorbar(cset1, ax=ax1, label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    ax1.legend(bbox_to_anchor=(-0.3, 1.05), loc='upper left', fontsize = legend_size)
    
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/dens_ecc_j_ang_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/dens_ecc_j_ang_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_obs_dens_vzp(parameters, variables, vzp, int_obs, int_j, int_k):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables

    fig, _axs = plt.subplots(nrows=1, ncols=3, figsize=(18, 5))
    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)
    axs = _axs.flatten()

    cset0 = axs[0].plot(int_k, int_j,     label="PEnGUIn")
    axs[0].plot(int_k,         int_obs, label="CI Tau Observed")
    cset1 = axs[1].pcolormesh(radius, theta, np.log10(density), vmin=-2, vmax=2)
    cset2 = axs[2].pcolormesh(radius, theta, vzp, vmin=vel_min, vmax=vel_max)

    axs[0].set_xlabel("Velocities [km/s]", fontsize = axis_size)
    axs[0].set_ylabel("Normalized Flux", fontsize = axis_size)
    axs[0].legend(fontsize = legend_size)
    axs[0].set_xlim(vel_min, vel_max)
    axs[0].set_ylim(lin_min, lin_max)
    axs[0].grid()

    axs[1].plot(radius[0,:], angle_1, color='k', linestyle='dashed', label="Inner Peri")
    axs[1].plot(radius[0,:], angle_2, color='r', linestyle='dashed', label="Outer Peri")
    #axs[1].plot(radius_1, theta[:,0], color='k', linestyle='dotted')
    #axs[1].plot(radius_2, theta[:,0], color='k', linestyle='dotted')

    axs[1].set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)
    axs[1].set_xlabel("Radius [au]", fontsize = axis_size)
    axs[1].set_xlim(RADIUS_IN, RADIUS_OUT)
    axs[1].set_ylabel("Theta [rads]", fontsize = axis_size)

    axs[2].set_title("PEnGUIn Inclined Velocity", fontsize = sub_title_size)
    axs[2].set_xlabel("Radius [au]", fontsize = axis_size)
    axs[2].set_xlim(RADIUS_IN, RADIUS_OUT)
    axs[2].set_ylabel("Theta [rads]", fontsize = axis_size)

    cbar1 = plt.colorbar(cset1, ax=axs[1], label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar1.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    axs[1].legend(bbox_to_anchor=(0.5, 0.95), loc='upper left', fontsize = legend_size)
    
    cbar2 = plt.colorbar(cset2, ax=axs[2], label="[km/s]")
    cbar2.set_label(label="[km/s]", fontsize = legend_size)
    
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        axs[0].set_title("Line Density", fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_vzp_" + str(k) + ".png")
    if run_type == "run_rotate":
        axs[0].set_title("Line Density", fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_vzp_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return
                     
                     
def plot_line_obs_dens_seg(parameters, variables, int_obs, int_j, int_k, int_b):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables
    
    b_len = len(int_b[:])
    b_label = []

    for i in range(b_len):
        b_label.append("Segment " + str(i))

    fig = plt.figure(figsize=(13, 5))
    ax0 = plt.subplot(121, projection='polar')
    ax1 = plt.subplot(122)

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)  

    ax0.grid(False)
    cset0 = ax0.pcolormesh(theta, radius, np.log10(density), vmin=-2, vmax=2)
    ax0.plot(angle_1, radius[0,:], color='k', linestyle='dashed', label="Inner Peri")
    ax0.plot(angle_2, radius[0,:], color='r', linestyle='dashed', label="Outer Peri")
    #ax0.plot(theta[:,0], radius_1, color = 'k', linestyle = 'dotted')
    #ax0.plot(theta[:,0], radius_2, color = 'k', linestyle = 'dotted')
    ax0.set_ylim(RADIUS_IN, RADIUS_OUT)
    ax0.set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)

    ax1.plot(int_k, int_j,        label="Sum of Segments")
    for i in range(b_len):
        ax1.plot(int_k, int_b[i], label=b_label[i])
    ax1.plot(int_k, int_obs,    label = "CI Tau Observed")
    ax1.set_xlabel("Velocities [km/s]", fontsize = axis_size)
    ax1.set_ylabel("Normalized Flux", fontsize = axis_size)
    ax1.set_xlim(vel_min, vel_max)
    ax1.set_ylim(lin_min, lin_max)
    ax1.grid()
    ax1.legend(fontsize = legend_size)

    cbar = plt.colorbar(cset0, ax=ax0, label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    ax0.legend(bbox_to_anchor=(-0.3, 1.05), loc='upper left', fontsize = legend_size)
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        ax1.set_title("Line Density", fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_seg_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax1.set_title("Line Density", fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_seg_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_dens_seg(parameters, variables, int_obs, int_j, int_k, int_b):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables
    
    b_len = len(int_b[:])
    b_label = []

    for i in range(b_len):
        b_label.append("Segment " + str(i))

    fig = plt.figure(figsize=(13, 5))
    ax0 = plt.subplot(121, projection='polar')
    ax1 = plt.subplot(122)

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)

    ax0.grid(False)
    cset0 = ax0.pcolormesh(theta, radius, np.log10(density), vmin=-2, vmax=2)
    ax0.plot(angle_1, radius[0,:], color='k', linestyle='dashed', label="Inner Peri")
    ax0.plot(angle_2, radius[0,:], color='r', linestyle='dashed', label="Outer Peri")
    #ax0.plot(theta[:,0], radius_1, color='k', linestyle='dotted')
    #ax0.plot(theta[:,0], radius_2, color='k', linestyle='dotted')
    ax0.set_ylim(RADIUS_IN, RADIUS_OUT)
    ax0.set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)
    
    ax1.plot(int_k, int_j,        label="Sum of Segments")
    for i in range(b_len):
        ax1.plot(int_k, int_b[i], label=b_label[i])
    ax1.set_xlabel("Velocities [km/s]", fontsize = axis_size)
    ax1.set_ylabel("Normalized Flux", fontsize = axis_size)
    ax1.set_xlim(vel_min, vel_max)
    ax1.set_ylim(lin_min, lin_max)
    ax1.grid()
    ax1.legend(fontsize=legend_size)

    cbar = plt.colorbar(cset0, ax=ax0, label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    ax0.legend(bbox_to_anchor=(-0.3, 1.05), loc='upper left', fontsize = legend_size)
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        ax1.set_title("Line Density", fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_dens_seg_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax1.set_title("Line Density", fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_dens_seg_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_obs_dens_wo_gap(parameters, variables, lines):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables
    int_obs, int_j, int_k, int_j_inner, int_j_gap, int_j_outer   = lines

    #int_j = int_j_inner + int_j_gap + int_j_outer
    max_int_j    = np.max(int_j)            #normalize calculated line to peak at 1
    int_j        = int_j       / max_int_j
    int_j_inner  = int_j_inner / max_int_j
    int_j_outer  = int_j_outer / max_int_j

    #Multi Plot - without Planet/Gap
    fig = plt.figure(figsize=(25, 10))
    ax0 = plt.subplot(121, projection='polar')
    ax1 = plt.subplot(122)

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)

    ax0.grid(False)
    cset0 = ax0.pcolormesh(theta, radius, np.log10(density), vmin=-2, vmax=2)
    ax0.plot(angle_1, radius[0,:], color='k', linestyle='dashed', label="Inner Peri")
    ax0.plot(angle_2, radius[0,:], color='r', linestyle='dashed', label="Outer Peri")
    #ax0.plot(theta[:,0], radius_1, color = 'k', linestyle = 'dotted')
    #ax0.plot(theta[:,0], radius_2, color = 'k', linestyle = 'dotted')
    ax0.set_ylim(RADIUS_IN, RADIUS_OUT)
    ax0.set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)

    ax1.plot(k_vel, int_j,        label = "Total")
    ax1.plot(k_vel, int_j_inner,  label = "Inner Disk")
    ax1.plot(k_vel, int_j_outer,  label = "Outer Disk")
    ax1.set_xlabel("Velocities [km/s]", fontsize = axis_size)
    ax1.set_ylabel("Normalized Flux", fontsize = axis_size)
    ax1.set_xlim(vel_min, vel_max)
    ax1.set_ylim(lin_min, lin_max)
    ax1.grid()
    ax1.legend(fontsize=legend_size)

    cbar = plt.colorbar(cset0, ax=ax0, label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    ax0.legend(bbox_to_anchor=(-0.3, 1.05), loc='upper left', fontsize = legend_size)
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        ax1.set_title("Line Density at orbit " + str(k), fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_dens_wo_gap_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax1.set_title("Line Density at orbit " + str(orbit_num) + " and frame " + str(frame), fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_dens_wo_gap_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_dens(parameters, variables, int_j):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables
    
    fig = plt.figure(figsize=(13, 5))
    ax0 = plt.subplot(121, projection='polar')
    ax1 = plt.subplot(122)

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)

    ax0.grid(False)
    cset0 = ax0.pcolormesh(theta, radius, np.log10(density), vmin=-2, vmax=2)
    ax0.plot(angle_1, radius[0,:], color='k', linestyle='dashed', label="Inner Peri")
    ax0.plot(angle_2, radius[0,:], color='r', linestyle='dashed', label="Outer Peri")
    #ax0.plot(theta[:,0], radius_1, color = 'k', linestyle = 'dotted')
    #ax0.plot(theta[:,0], radius_2, color = 'k', linestyle = 'dotted')
    ax0.set_ylim(RADIUS_IN, RADIUS_OUT)
    ax0.set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)
    
    
    ax1.plot(k_vel, int_j,        label = "PEnGUIn")
    ax1.set_xlabel("Velocities [km/s]", fontsize = axis_size)
    ax1.set_ylabel("Normalized Flux", fontsize = axis_size)
    ax1.set_xlim(vel_min, vel_max)
    ax1.set_ylim(lin_min, lin_max)
    ax1.grid()
    ax1.legend(fontsize=legend_size)

    cbar = plt.colorbar(cset0, ax=ax0, label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    ax0.legend(bbox_to_anchor=(-0.3, 1.05), loc='upper left', fontsize = legend_size)
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        ax1.set_title("Line Density at orbit " + str(k), fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_dens_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax1.set_title("Line Density at orbit " + str(orbit_num) + " and frame " + str(frame), fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_dens_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_obs_dens(parameters, variables, int_obs, int_j, int_k):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables

    fig = plt.figure(figsize=(13, 5))
    ax0 = plt.subplot(121, projection='polar')
    ax1 = plt.subplot(122)

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)

    ax0.grid(False)
    cset0 = ax0.pcolormesh(theta, radius, np.log10(density), vmin=-2, vmax=2)
    ax0.plot(angle_1, radius[0,:], color='k', linestyle='dashed', label="Inner Peri")
    ax0.plot(angle_2, radius[0,:], color='r', linestyle='dashed', label="Outer Peri")
    #ax0.plot(theta[:,0], radius_1, color = 'k', linestyle = 'dotted')
    #ax0.plot(theta[:,0], radius_2, color = 'k', linestyle = 'dotted')
    ax0.set_ylim(0,0.3)
    ax0.set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)
    
    ax1.plot(int_k,int_obs, color='orange', label="CI Tau Observed")
    ax1.plot(int_k,int_j, color='b', label="PEnGUIUn")
    ax1.set_xlabel("Velocities [km/s]", fontsize = axis_size)
    ax1.set_ylabel("Normalized Flux", fontsize = axis_size)
    ax1.set_xlim(vel_min, vel_max)
    ax1.set_ylim(lin_min, lin_max)
    ax1.grid()
    ax1.legend(fontsize=legend_size)

    cbar = plt.colorbar(cset0, ax=ax0, label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    ax0.legend(bbox_to_anchor=(-0.3, 1.05), loc='upper left', fontsize = legend_size)
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        ax1.set_title("Density at orbit " + str(k), fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax1.set_title("Line Density at orbit " + str(orbit_num) + " and frame " + str(frame), fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_obs_std_dens(parameters, variables, int_obs, int_j, int_k, norm_std):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, radius_1, radius_2, angle_1, angle_2 = variables

    std1_upper = int_obs + norm_std
    std1_lower = int_obs - norm_std
    std2_upper = int_obs + 2*norm_std
    std2_lower = int_obs - 2*norm_std

    fig = plt.figure(figsize=(13, 5))
    ax0 = plt.subplot(121, projection='polar')
    ax1 = plt.subplot(122)

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)

    ax0.grid(False)
    cset0 = ax0.pcolormesh(theta, radius, np.log10(density), vmin=-2, vmax=2)
    ax0.plot(angle_1, radius[0,:], color='k', linestyle='dashed', label="Inner Peri")
    ax0.plot(angle_2, radius[0,:], color='r', linestyle='dashed', label="Outer Peri")
    #ax0.plot(theta[:,0], radius_1, color = 'k', linestyle = 'dotted')
    #ax0.plot(theta[:,0], radius_2, color = 'k', linestyle = 'dotted')
    ax0.set_ylim(RADIUS_IN, RADIUS_OUT)
    ax0.set_title("Calculated PEnGUIn Densities", fontsize = sub_title_size)
    
    ax1.plot(int_k, int_obs, color='orange', label="CI Tau Observed")
    ax1.fill_between(int_k, std1_lower, std1_upper, color='g', alpha=0.6, label="1 STD")
    ax1.fill_between(int_k, std2_lower, std2_upper, color='g', alpha=0.3, label="2 STD")
    ax1.plot(int_k, int_j, color='b', label="PEnGUIUn")
    ax1.set_xlabel("Velocities [km/s]", fontsize = axis_size)
    ax1.set_ylabel("Normalized Flux", fontsize = axis_size)
    ax1.set_xlim(vel_min, vel_max)
    ax1.set_ylim(lin_min, lin_max)
    ax1.grid()
    ax1.legend(fontsize=legend_size)

    cbar = plt.colorbar(cset0, ax=ax0, label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$')
    cbar.set_label(label=r'$\log_{10}\left(\frac{\Sigma}{\Sigma_0}\right)$', fontsize = legend_size)
    ax0.legend(bbox_to_anchor=(-0.3, 1.05), loc='upper left', fontsize = legend_size)
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        ax1.set_title("Line Density", fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_std_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax1.set_title("Line Density at orbit " + str(orbit_num) + " and frame " + str(frame), fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_std_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


#-------------------------------------------------------------------------------------------------------
#other plots -------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def plot_col_dens(parameters, varibles, params2):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, density, temp, part_func     = varibles
    vel_j, ein_j, eng_j, g1j, g0j        = params2
    
    removed = np.transpose(density) * DENS_FACTOR
    removed[400:500,1050:1110] = removed[400:500,1050:1110] - removed[400:500,1050:1110]
    
    constant1 = 1 / (8*B*100000*(np.pi)**(3./2.))

    dens = np.average(removed,axis=1)
    N0_detected  = N0
    N0_simulated = N0 * 100.0
    N_detected  = N0_detected  * (radius/RADIUS_IN)**BETA
    N_simulated = N0_simulated * (radius/RADIUS_IN)**BETA
    N_factor    = DENS_FACTOR  * (radius/RADIUS_IN)**BETA
    thick_dense_P42 = part_func /( np.exp(-eng_j[0]  / (BOLTZ * temp)) * (ein_j[0]  * g0j[0]**2  / (g1j[0]  * vel_j[0]**3))  * constant1)
    thick_dense_P20 = part_func /( np.exp(-eng_j[-1] / (BOLTZ * temp)) * (ein_j[-1] * g0j[-1]**2  / (g1j[-1] * vel_j[-1]**3)) * constant1)


    fig, ax =plt.subplots()
    ax.plot(radius, np.log10(dens),                 color='b', label="AVG PEnGUIn")
    ax.plot(radius, np.log10(N_detected),           color="g", label="Detected")
    ax.plot(radius, np.log10(N_simulated),          color="orange", label="Simulated")
    ax.plot(radius, np.log10(thick_dense_P42[0,:]), color="r",  linestyle="dotted", label="Optically Thick P42 Sur")
    ax.plot(radius, np.log10(thick_dense_P20[0,:]), color="r",  linestyle="dashed", label="Optically Thick P20 Sur")
    plt.legend(fontsize = legend_size)
    plt.title("Column Density at Orbit " +str(k), fontsize = title_size)
    plt.xlabel("Radius [au]", fontsize = axis_size)
    plt.ylabel("Column Density [log(cm-2)]", fontsize = axis_size)
    plt.xlim(RADIUS_IN, 0.35)
    plt.ylim(16,24)
    if run_type == "run_stationary":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Col Dens" + "/col_dens_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Col Dens" + "/col_dens_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return

    
def plot_opt_depth(parameters, varibles):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, optical_t = varibles
    
    removed = np.transpose(optical_t)
    removed[400:500,1050:1110] = removed[400:500,1050:1110] - removed[400:500,1050:1110]
    
    optical_depth = np.average(removed, axis=1)

    fig, ax = plt.subplots()
    ax.plot(radius, optical_depth, label="PEN Optical")
    plt.legend(fontsize = legend_size)
    plt.xlabel("Radius [au]", fontsize = axis_size)
    plt.ylabel("Optical Depth", fontsize = axis_size)
    plt.title("Optical Depth Average at orbit " + str(orbit_num), fontsize = title_size)
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Optical Depth" + "/optical_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Optical Depth" + "/optical_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_opt_depth_all(parameters, varibles):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, optical_t = varibles
    
    removed = np.transpose(optical_t)
    removed[:,400:500,1050:1110] = removed[:,400:500,1050:1110] - removed[:,400:500,1050:1110]

    optical_depth = np.average(removed, axis=1)

    fig, ax = plt.subplots()
    for i in range(len(optical_depth[0,:])):
        ax.plot(radius[0,:], optical_depth[:,i], label=i)
    plt.xlabel("Radius [au]", fontsize = axis_size)
    plt.ylabel("Optical Depth", fontsize = axis_size)
    plt.title("Optical Depth Average at orbit " + str(orbit_num), fontsize = title_size)
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Optical Depth" + "/optical_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Optical Depth" + "/optical_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return



def plot_line_vzp(parameters, varibles):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, theta, density, int_j, radius_1, radius_2     = varibles
    

    #Multi Plot - polar density and line
    fig = plt.figure(figsize=(13, 5))
    ax0 = plt.subplot(121, projection='polar')
    ax1 = plt.subplot(122)

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k), fontsize = title_size)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' and at Frame ' + str(frame), fontsize = title_size)

    ax0.grid(False)
    cset0 = ax0.pcolormesh(theta, radius, vzp, vmin=120, vmax=130)
    ax0.set_ylim(0,0.3)
    ax0.set_title("Calculated INclined Velcoties", fontsize = sub_title_size)
    
    #ax0.plot(theta[:,0], radius_1, color = 'k', linestyle = 'dotted')
    #ax0.plot(theta[:,0], radius_2, color = 'k', linestyle = 'dotted')

    ax1.plot(k_vel, int_j,        label = "PEnGUIn")
    ax1.set_xlabel("Velocities [km/s]", fontsize = axis_size)
    ax1.set_ylabel("Normalized Flux", fontsize = axis_size)
    ax1.set_xlim(-200, 200)
    ax1.grid()
    ax1.legend(fontsize=legend_size)

    fig.colorbar(cset0, ax=ax0, label="[log(cm-2)]")
    plt.tight_layout(pad=padding)
    
    if run_type == "run_stationary":
        ax1.set_title("Line Density at orbit " + str(k), fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_vzp_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax1.set_title("Line Density at orbit " + str(orbit_num) + " and frame " + str(frame), fontsize = sub_title_size)
        plt.savefig("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_vzp_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return