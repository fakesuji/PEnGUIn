"""
Author:  Cory Padgett
Advisor: Dr. Jeffrey Fung
Email:   cpadge4@clemson.edu
"""

#imports
import numpy as np
import matplotlib.pyplot as plt


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
MASS_STAR   = 1.02 * MASS_SOL   #g
MASS_PLANET = 11.6 * MASS_JUP   #g
INCLIN  = np.radians(71.0)      #radians
T0      = 2770           #K
ALPHA   = -0.33          #
N0      = 10**19         #
BETA    = -2.0           #
R_IN    = 0.048          #au
B       = 8.8            #km/s
B_CGS   = B * 100000     #cm/s


#-------------------------------------------------------------------------------------------------------
#plotting functions
#-------------------------------------------------------------------------------------------------------
#plot line density spectrum ----------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def plot_line(parameters, k_vel, int_j):
    #plot ro-vibrational line density
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters

    max_int_j = np.max(int_j)      #normalize calculated line to peak at 1
    int_j     = int_j / max_int_j

    fig, ax = plt.subplots()
    ax.plot(k_vel,int_j, label="PEnGUIUn")
    plt.legend()
    plt.xlabel("Velocity [km/s]")
    plt.ylabel("Normalized Flux")
    if run_type == "run_stationary":
        plt.title("Stationary Line Density at orbit " + str(k))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.title("Rotational Line Density at orbit " + str(orbit_num) + " and frame " + str(frame))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_obs_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_obs(parameters, k_vel, int_j):
    #plot ro-vibrational line density with observation line overlay
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters

    int_j_obs = np.loadtxt("Data/CITau_obs_line.txt") / 0.95   #0.95 factor to bring observation peaks higher adjusting for noise
    max_int_j = np.max(int_j)      #normalize calculated line to peak at 1
    int_j     = int_j / max_int_j

    fig, ax = plt.subplots()
    ax.plot(k_vel,int_j,     label="PEnGUIUn")
    ax.plot(k_vel,int_j_obs, label="CI Tau Observed")
    plt.legend()
    plt.xlabel("Velocity [km/s]")
    plt.ylabel("Normalized Flux")
    if run_type == "run_stationary":
        plt.title("Stationary Line Density at orbit " + str(k))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_obs_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.title("Rotational Line Density at orbit " + str(orbit_num) + " and frame " + str(frame))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_obs_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_seg(parameters, varibles):
    #plot ro-vibrational line density and in spatial segments
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    int_j_inner, int_j_gap, int_j_outer, k_vel = varibles

    int_j = int_j_inner + int_j_gap + int_j_outer
    max_int_j    = np.max(int_j)           #normalize calculated line to peak at 1
    int_j        = int_j / max_int_j
    int_j_inner  = int_j_inner / max_int_j
    int_j_outer  = int_j_outer / max_int_j
    int_j_gap    = int_j_gap / max_int_j

    fig, ax =plt.subplots()
    ax.plot(k_vel,int_j_inner,  label="PEnGUIUn  Inner Disk")
    ax.plot(k_vel,int_j_gap,    label="PEnGUIUn  Gap")
    ax.plot(k_vel,int_j_outer,  label="PEnGUIUn  Outer Disk")
    ax.plot(k_vel,int_j,        label="PEnGUIUn  Sum")
    plt.legend()
    plt.xlabel("Velocity [km/s]")
    plt.ylabel("Normalized Flux")
    if run_type == "run_stationary":
        plt.title("Stationary Line Density at orbit " + str(k))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_obs_segs_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.title("Rotational Line Density at orbit " + str(orbit_num) + " and frame " + str(frame))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_obs_segs" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_obs_seg(parameters, varibles):
    #plot ro-vibrational line density with observation line overlay and in spatial segments
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    int_j_inner, int_j_gap, int_j_outer, k_vel = varibles

    int_j_obs = np.loadtxt("Data/CITau_obs_line.txt") / 0.95   #0.95 factor to bring observation peaks higher adjusting for noise
    int_j = int_j_inner + int_j_gap + int_j_outer
    max_int_j    = np.max(int_j)           #normalize calculated line to peak at 1
    int_j        = int_j / max_int_j
    int_j_inner  = int_j_inner / max_int_j
    int_j_outer  = int_j_outer / max_int_j
    int_j_gap    = int_j_gap / max_int_j

    fig, ax =plt.subplots()
    ax.plot(k_vel,int_j_inner,  label="PEnGUIUn Inner Disk")
    ax.plot(k_vel,int_j_gap,    label="PEnGUIUn Gap")
    ax.plot(k_vel,int_j_outer,  label="PEnGUIUn Outer Disk")
    ax.plot(k_vel,int_j,        label="PEnGUIUn Total")
    ax.plot(k_vel,int_j_obs,    label="CI Tau Observed")
    plt.legend()
    plt.xlabel("Velocity [km/s]")
    plt.ylabel("Normalized Flux")
    if run_type == "run_stationary":
        plt.title("Stationary Line Density at orbit " + str(k))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_obs_segs_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.title("Rotational Line Density at orbit " + str(orbit_num) + " and frame " + str(frame))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Line Density" + "/line_obs_segs" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


#-------------------------------------------------------------------------------------------------------
#polar plots -------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def plot_vel(parameters, radius, theta, vel, vel_type):
    #Polar velocity plot
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6,5))
    ax.grid(False)
    plt.pcolormesh(theta, radius, vel)
    ax.set_rmax(0.3)
    ax.set_rticks([0.1, 0.2, 0.3])  # Less radial ticks
    ax.set_rlabel_position(-22.5)   # Move radial labels away from plotted line
    plt.colorbar(label="[km/s]")
    if run_type == "run_stationary":
        ax.set_title("PEnGUIn " + vel_type + " Velocities at time " + str(time) + " and Orbit " + str(k))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + vel_type + "/vel_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax.set_title("PEnGUIn " + vel_type + " Velocities at time " + str(time) + ",\n Orbit " + str(orbit_num) + ", and Frame " + str(frame))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + vel_type + "/vel_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_pol_dens(parameters, radius, theta, radius_1, radius_2, density):
    #Polar Density Map
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6,5))
    ax.grid(False)
    plt.pcolormesh(theta, radius, np.log10(density), vmin=15, vmax=20)
    plt.plot(theta[:,0], radius_1, color='k', linestyle='dotted')     #plot visual inner gap boundary
    plt.plot(theta[:,0], radius_2, color='k', linestyle='dotted')     #plot visual outer gap boundary
    ax.set_rmax(0.3)
    ax.set_rticks([0.1, 0.2, 0.3])  # Less radial ticks
    ax.set_rlabel_position(-22.5)   # Move radial labels away from plotted line
    plt.colorbar(label="[log(cm-2)]")
    if run_type == "run_stationary":
        ax.set_title("Calculated PEnGUIn Densities at time " + str(time) + " and Orbit " + str(k))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Polar Density" + "/pol_dens_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax.set_title("Calculated PEnGUIn Densities at time " + str(time) + ",\n Orbit " + str(orbit_num) + ", and Frame " + str(frame))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Polar Density" + "/pol_dens_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_pol_vzp(parameters, radius, theta, radius_1, radius_2, vzp):
    #Polar Inclined Velocity Map
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.grid(False)
    plt.pcolormesh(theta, radius, vzp, vmin=-150.0, vmax=150.0)
    plt.plot(theta[:,0], radius_1, color='k', linestyle='dotted')     #plot visual inner gap boundary
    plt.plot(theta[:,0], radius_2, color='k', linestyle='dotted')     #plot visual outer gap boundary
    ax.set_rmax(0.3)
    ax.set_rticks([0.1, 0.2, 0.3])  # Less radial ticks
    ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    plt.colorbar(label="[km/s]")
    if run_type == "run_stationary":
        ax.set_title("Projected PEnGUIn Velocities at time " + str(time) + " and Orbit " + str(k))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Proj-Pol Velocity" + "/pol_vzp_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax.set_title("Projected PEnGUIn Velocities at time " + str(time) + ",\n Orbit " + str(orbit_num) + ", and Frame " + str(frame))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Proj-Pol Velocity" + "/pol_vzp_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_pol_temp(parameters, radius, theta, radius_1, radius_2, temp):
    #Polar Temperature Map
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.grid(False)
    plt.pcolormesh(theta, radius, np.log10(temp))
    plt.plot(theta[:,0], radius_1, color='k', linestyle='dotted')     #plot visual inner gap boundary
    plt.plot(theta[:,0], radius_2, color='k', linestyle='dotted')     #plot visual outer gap boundary
    ax.set_rmax(0.3)
    ax.set_rticks([0.1, 0.2, 0.3])  # Less radial ticks
    ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    plt.colorbar(label="[log10(k)]")
    if run_type == "run_stationary":
        ax.set_title("Surface Temperatures at time " + str(time) + " and Orbit " + str(k))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Temperatures" + "/temp_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax.set_title("Surface Temperatures at time " + str(time) + ",\n Orbit " + str(orbit_num) + ", and Frame " + str(frame))
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Temperatures" + "/temp_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


#-------------------------------------------------------------------------------------------------------
#Multi-info Plots --------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def plot_line_obs_dens_vzp(parameters, varibles):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    k_vel, int_j, radius, theta, density, vzp, radius_1, radius_2, R_IN = varibles

    int_j_obs = np.loadtxt("Data/CITau_obs_line.txt") / 0.95   #0.95 factor to bring observation peaks higher adjusting for noise
    max_int_j = np.max(int_j)       #normalize calculated line to peak at 1
    int_j     = int_j / max_int_j

    fig, _axs = plt.subplots(nrows=1, ncols=3, figsize=(18, 5))
    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k) + ' at Time ' + str(time))
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' at Frame ' + str(frame) + ' at Time ' + str(time))
    axs = _axs.flatten()

    cset0 = axs[0].plot(k_vel, int_j,     label="PEnGUIn")
    axs[0].plot(k_vel,         int_j_obs, label="CI Tau Observed")
    cset1 = axs[1].pcolormesh(radius, theta, np.log10(density), vmin=15, vmax=20)
    cset2 = axs[2].pcolormesh(radius, theta, vzp, vmin=-150, vmax=150)

    axs[0].set_xlabel("Velocities [km/s]")
    axs[0].set_ylabel("Normalized Flux")
    axs[0].legend()

    axs[1].set_title("Calculated PEnGUIn Densities")
    axs[1].set_xlabel("Radius [au]")
    axs[1].set_xlim(R_IN, 0.5)
    axs[1].set_ylabel("Theta [rads]")

    axs[2].set_title("PEnGUIn Inclined Velocity")
    axs[2].set_xlabel("Radius [au]")
    axs[2].set_xlim(R_IN, 0.5)
    axs[2].set_ylabel("Theta [rads]")

    axs[1].plot(radius_1, theta[:,0], color='k', linestyle='dotted')     #plot visual inner gap boundary
    axs[1].plot(radius_2, theta[:,0], color='k', linestyle='dotted')     #plot visual outer gap boundary

    fig.colorbar(cset1, ax=axs[1], label="[log(cm-2)]")
    fig.colorbar(cset2, ax=axs[2], label="[km/s]")
    if run_type == "run_stationary":
        axs[0].set_title("Stationary Line Density at orbit " + str(k), fontsize=15)
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_vzp_" + str(k) + ".png")
    if run_type == "run_rotate":
        axs[0].set_title("Rotational Line Density at orbit " + str(orbit_num) + " and frame " + str(frame), fontsize=15)
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_vzp_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return
                     
                     
def plot_line_obs_dens_seg(parameters, varibles, lines):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    k_vel, radius, theta, density, radius_1, radius_2    = varibles
    int_j_inner, int_j_gap, int_j_outer                  = lines

    int_j_obs = np.loadtxt("Data/CITau_obs_line.txt") / 0.95   #0.95 factor to bring observation peaks higher adjusting for noise
    int_j = int_j_inner + int_j_gap + int_j_outer
    max_int_j    = np.max(int_j)            #normalize calculated line to peak at 1
    int_j        = int_j       / max_int_j
    int_j_inner  = int_j_inner / max_int_j
    int_j_outer  = int_j_outer / max_int_j
    int_j_gap    = int_j_gap   / max_int_j

    #Multi Plot
    fig = plt.figure(figsize=(25, 10))
    ax0 = plt.subplot(121, projection='polar')
    ax1 = plt.subplot(122)

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k) + ' at Time ' + str(time), fontsize=16)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' at Frame ' + str(frame) + ' at Time ' + str(time), fontsize=16)    #k-frame_start

    ax0.grid(False)
    cset0 = ax0.pcolormesh(theta, radius, np.log10(np.transpose(density)), vmin=15, vmax=20)
    ax0.set_ylim(0,0.3)
    ax0.set_title("Calculated PEnGUIn Densities", fontsize=15)
    ax0.plot(theta[:,0], radius_1, color = 'k', linestyle = 'dotted')     #plot visual inner gap boundary
    ax0.plot(theta[:,0], radius_2, color = 'k', linestyle = 'dotted')     #plot visual outer gap boundary

    ax1.plot(k_vel, int_j,        label = "PEN Inc")
    ax1.plot(k_vel, int_j_inner,  label = "Inner Disk")
    ax1.plot(k_vel, int_j_gap,    label = "Gap w Planet")
    ax1.plot(k_vel, int_j_outer,  label = "Outer Disk")
    ax1.plot(k_vel, int_j_obs,    label = "CI Tau Observed")
    ax1.set_xlabel("Velocities [km/s]", fontsize=12)
    ax1.set_ylabel("Normalized Flux",   fontsize=12)
    ax1.set_xlim(-200, 200)
    ax1.legend(fontsize = 15)

    fig.colorbar(cset0, ax=ax0, label="[log(cm-2)]")
    if run_type == "run_stationary":
        ax1.set_title("Stationary Line Density at orbit " + str(k), fontsize=15)
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax1.set_title("Rotational Line Density at orbit " + str(orbit_num) + " and frame " + str(frame), fontsize=15)
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_dens_seg(parameters, varibles, lines):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    k_vel, radius, theta, density, radius_1, radius_2    = varibles
    int_j_inner, int_j_gap, int_j_outer                  = lines

    int_j = int_j_inner + int_j_gap + int_j_outer
    max_int_j    = np.max(int_j)            #normalize calculated line to peak at 1
    int_j        = int_j       / max_int_j
    int_j_inner  = int_j_inner / max_int_j
    int_j_outer  = int_j_outer / max_int_j
    int_j_gap    = int_j_gap   / max_int_j

    #Multi Plot
    fig = plt.figure(figsize=(25, 10))
    ax0 = plt.subplot(121, projection='polar')
    ax1 = plt.subplot(122)

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k) + ' at Time ' + str(time), fontsize=16)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' at Frame ' + str(frame) + ' at Time ' + str(time), fontsize=16)

    ax0.grid(False)
    cset0 = ax0.pcolormesh(theta, radius, np.log10(np.transpose(density)), vmin=15, vmax=20)
    ax0.set_ylim(0,0.3)
    ax0.set_title("Calculated PEnGUIn Densities", fontsize=15)
    ax0.plot(theta[:,0], radius_1, color='k', linestyle='dotted')     #plot visual inner gap boundary
    ax0.plot(theta[:,0], radius_2, color='k', linestyle='dotted')     #plot visual outer gap boundary

    ax1.plot(k_vel, int_j,        label="Total")
    ax1.plot(k_vel, int_j_inner,  label="Inner Disk")
    ax1.plot(k_vel, int_j_gap,    label="Gap w Planet")
    ax1.plot(k_vel, int_j_outer,  label="Outer Disk")
    ax1.set_xlabel("Velocities [km/s]", fontsize=12)
    ax1.set_ylabel("Normalized Flux",   fontsize=12)
    ax1.set_xlim(-200,200)
    ax1.legend(fontsize=15)

    fig.colorbar(cset0, ax=ax0, label="[log(cm-2)]")
    if run_type == "run_stationary":
        ax1.set_title("Stationary Line Density at orbit " + str(k), fontsize=15)
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_dens_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax1.set_title("Rotational Line Density at orbit " + str(orbit_num) + " and frame " + str(frame), fontsize=15)
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_dens_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_obs_dens_wo_gap(parameters, varibles, lines):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    k_vel, radius, theta, density, radius_1, radius_2    = varibles
    int_j_inner, int_j_gap, int_j_outer                  = lines

    int_j_obs = np.loadtxt("Data/CITau_obs_line.txt") / 0.95   #0.95 factor to bring observation peaks higher adjusting for noise
    int_j = int_j_inner + int_j_gap + int_j_outer
    max_int_j    = np.max(int_j)            #normalize calculated line to peak at 1
    int_j        = int_j       / max_int_j
    int_j_inner  = int_j_inner / max_int_j
    int_j_outer  = int_j_outer / max_int_j

    #Multi Plot - without Planet/Gap
    fig = plt.figure(figsize=(25, 10))
    ax0 = plt.subplot(121, projection='polar')
    ax1 = plt.subplot(122)

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k) + ' at Time ' + str(time), fontsize=16)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' at Frame ' + str(frame) + ' at Time ' + str(time), fontsize=16)

    ax0.grid(False)
    cset0 = ax0.pcolormesh(theta, radius, np.log10(density), vmin=15, vmax=20)
    ax0.set_ylim(0,0.3)
    ax0.set_title("Calculated PEnGUIn Densities", fontsize=15)
    ax0.plot(theta[:,0], radius_1, color = 'k', linestyle = 'dotted')     #plot visual inner gap boundary
    ax0.plot(theta[:,0], radius_2, color = 'k', linestyle = 'dotted')     #plot visual outer gap boundary

    ax1.plot(k_vel, int_j,        label = "Total")
    ax1.plot(k_vel, int_j_inner,  label = "Inner Disk")
    ax1.plot(k_vel, int_j_outer,  label = "Outer Disk")
    ax1.set_xlabel("Velocities [km/s]", fontsize=12)
    ax1.set_ylabel("Normalized Flux",   fontsize=12)
    ax1.set_xlim(-200,200)
    ax1.legend(fontsize=15)

    fig.colorbar(cset0, ax=ax0, label="[log(cm-2)]")
    if run_type == "run_stationary":
        ax1.set_title("Stationary Line Density at orbit " + str(k), fontsize=15)
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_dens_wo_gap_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax1.set_title("Rotational Line Density at orbit " + str(orbit_num) + " and frame " + str(frame), fontsize=15)
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_dens_wo_gap_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_dens(parameters, varibles):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    k_vel, radius, theta, density, int_j, radius_1, radius_2     = varibles
    
    max_int_j = np.max(int_j)       #normalize calculated line to peak at 1
    int_j     = int_j / max_int_j

    #Multi Plot - polar density and line
    fig = plt.figure(figsize=(12, 5))
    ax0 = plt.subplot(121, projection='polar')
    ax1 = plt.subplot(122)

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k) + ' at Time ' + str(time), fontsize=16)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' at Frame ' + str(frame) + ' at Time ' + str(time), fontsize=16)

    ax0.grid(False)
    cset0 = ax0.pcolormesh(theta, radius, np.log10(density), vmin=15, vmax=20)
    ax0.set_ylim(0,0.3)
    ax0.set_title("Calculated PEnGUIn Densities", fontsize=15)
    ax0.plot(theta[:,0], radius_1, color = 'k', linestyle = 'dotted')     #plot visual inner gap boundary
    ax0.plot(theta[:,0], radius_2, color = 'k', linestyle = 'dotted')     #plot visual outer gap boundary

    ax1.plot(k_vel, int_j,        label = "PEnGUIn")
    ax1.set_xlabel("Velocities [km/s]", fontsize=12)
    ax1.set_ylabel("Normalized Flux",   fontsize=12)
    ax1.set_xlim(-200,200)
    ax1.legend(fontsize=12)

    fig.colorbar(cset0, ax=ax0, label="[log(cm-2)]")
    if run_type == "run_stationary":
        ax1.set_title("Stationary Line Density at orbit " + str(k), fontsize=15)
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_dens_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax1.set_title("Rotational Line Density at orbit " + str(orbit_num) + " and frame " + str(frame), fontsize=15)
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_dens_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


def plot_line_obs_dens(parameters, varibles):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    k_vel, radius, theta, density, int_j, radius_1, radius_2     = varibles

    int_j_obs = np.loadtxt("Data/CITau_obs_line.txt") / 0.95   #0.95 factor to bring observation peaks higher adjusting for noise
    max_int_j = np.max(int_j)       #normalize calculated line to peak at 1
    int_j     = int_j / max_int_j

    #Multi Plot - polar density and line
    fig = plt.figure(figsize=(12, 5))
    ax0 = plt.subplot(121, projection='polar')
    ax1 = plt.subplot(122)

    if run_type == "run_stationary":
        fig.suptitle('Orbit ' +str(k) + ' at Time ' + str(time), fontsize=16)
    if run_type == "run_rotate":
        fig.suptitle('Orbit ' +str(orbit_num) + ' at Frame ' + str(frame) + ' at Time ' + str(time), fontsize=16)    #k-frame_start

    ax0.grid(False)
    cset0 = ax0.pcolormesh(theta, radius, np.log10(density), vmin=15, vmax=20)
    ax0.set_ylim(0,0.3)
    ax0.set_title("Calculated PEnGUIn Densities", fontsize=15)
    ax0.plot(theta[:,0], radius_1, color = 'k', linestyle = 'dotted')     #plot visual inner gap boundary
    ax0.plot(theta[:,0], radius_2, color = 'k', linestyle = 'dotted')     #plot visual outer gap boundary

    ax1.plot(k_vel, int_j,        label = "PEnGUIn")
    ax1.plot(k_vel, int_j_obs,    label = "Observation")
    ax1.set_xlabel("Velocities [km/s]", fontsize=12)
    ax1.set_ylabel("Normalized Flux",   fontsize=12)
    ax1.set_xlim(-200,200)
    ax1.legend(fontsize=10)

    fig.colorbar(cset0, ax=ax0, label="[log(cm-2)]")
    if run_type == "run_stationary":
        ax1.set_title("Stationary Line Density at orbit " + str(k), fontsize=15)
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_" + str(k) + ".png")
    if run_type == "run_rotate":
        ax1.set_title("Rotational Line Density at orbit " + str(orbit_num) + " and frame " + str(frame), fontsize=15)
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Multi Plot" + "/line_obs_dens_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return


#-------------------------------------------------------------------------------------------------------
#other plots -------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def plot_col_dens(parameters, varibles, params2):
    run_type, run_name, analysis_name, time, k, frame, orbit_num = parameters
    radius, density, temp, part_func, dens_factor     = varibles
    vel_j, ein_j, eng_j, g1j, g0j                     = params2
    
    removed = np.transpose(density)
    removed[400:500,1050:1110] = removed[400:500,1050:1110] - removed[400:500,1050:1110]
    
    constant1 = 1 / (8*B*100000*(np.pi)**(3./2.))

    dens = np.average(removed,axis=1)
    N0_detected  = N0
    N0_simulated = N0 * 100.0
    N_detected  = N0_detected  * (radius/R_IN)**BETA
    N_simulated = N0_simulated * (radius/R_IN)**BETA
    N_factor    = dens_factor  * (radius/R_IN)**BETA
    thick_dense_P42 = part_func /( np.exp(-eng_j[0]  / (BOLTZ * temp)) * (ein_j[0]  * g0j[0]**2  / (g1j[0]  * vel_j[0]**3))  * constant1)
    thick_dense_P20 = part_func /( np.exp(-eng_j[-1] / (BOLTZ * temp)) * (ein_j[-1] * g0j[-1]**2  / (g1j[-1] * vel_j[-1]**3)) * constant1)


    fig, ax =plt.subplots()
    ax.plot(radius, np.log10(dens),                 color='b', label="AVG PEnGUIn")
    ax.plot(radius, np.log10(N_detected),      color="g", label="detected")
    ax.plot(radius, np.log10(N_simulated),     color="orange", label="simulated")
    ax.plot(radius, np.log10(thick_dense_P42[0,:]), color="r", linestyle="dotted", label="Optically Thick P42 Sur")
    ax.plot(radius, np.log10(thick_dense_P20[0,:]), color="r", linestyle="dashed", label="Optically Thick P20 Sur")
    plt.legend()
    plt.title("Column Density at time " + str(time) + " and orbit " +str(k))
    plt.xlabel("Radius [au]")
    plt.ylabel("Column Density [log(cm-2)]")
    plt.xlim(R_IN, 0.35)
    plt.ylim(16,24)
    if run_type == "run_stationary":
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Col Dens" + "/col_dens_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Col Dens" + "/col_dens_" + str(orbit_num) + "_" + str(frame) + ".png")
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
    plt.legend()
    plt.xlabel("Radius [au]")
    plt.ylabel("Optical Depth")
    plt.title("Optical Depth Average at time " + str(time))
    if run_type == "run_stationary":
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Optical Depth" + "/optical_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Optical Depth" + "/optical_" + str(orbit_num) + "_" + str(frame) + ".png")
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
    plt.xlabel("Radius [au]")
    plt.ylabel("Optical Depth")
    plt.title("Optical Depth Average at time " + str(time))
    if run_type == "run_stationary":
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Optical Depth" + "/optical_" + str(k) + ".png")
    if run_type == "run_rotate":
        plt.savefig("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + "Optical Depth" + "/optical_" + str(orbit_num) + "_" + str(frame) + ".png")
    plt.close()
    return