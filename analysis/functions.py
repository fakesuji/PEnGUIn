"""
Author:  Cory Padgett
Advisor: Dr. Jeffrey Fung
Email:   cpadge4@clemson.edu
"""

#imports
import os
from datetime import datetime
import csv
import numpy as np
import read_penguin as r
from scipy.stats import circmean, circvar, circstd
from scipy.stats import linregress

#-------------------------------------------------------------------------------------------------------
#user info
#-------------------------------------------------------------------------------------------------------
LOAD_PATH  = "/scratch/cpadge4/"

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
H_MASS     = 2.01568 #g/mol
ECC_NORM   = 1/np.pi


#system constants
#---------------------------------------
#planet-star constants
MASS_STAR      = 1.02 * MASS_SOL   #g
MASS_PLANET    = 2.0  * MASS_JUP   #g
PLANET_SCALING = 0.164
ORBIT_FREQ     = 15.0568536075
orbit_interval = 10.0
#---------------------------------------
#disk constants
INCLIN  = np.radians(71.0)      #radians
T0      = 2770           #K
ALPHA   = -0.33          #
N0      = 2.0E19         #
BETA    = -2.0           #
B       = 8.8            #km/s
B_CGS   = B * 100000     #cm/s
INNER   = 0.8 * PLANET_SCALING
OUTER   = 1.2 * PLANET_SCALING
RADIUS_IN    =  0.25 * PLANET_SCALING
RADIUS_OUT   =  10.0 * PLANET_SCALING
#---------------------------------------


#arrays
k_vel = np.linspace(-300,300,601)


#-------------------------------------------------------------------------------------------------------
#functions
#-------------------------------------------------------------------------------------------------------
#loading data ------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def read_csv(file_path):
    with open(file_path, 'r') as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader)
        data = list(zip(*[map(float, row) for row in reader]))

    column_data = {header[i]: np.array(data[i]) for i in range(len(header))}
    return column_data

def get_stats_data(RUN_NAME, ANALYSIS_NAME):
    file_name = ["IJ_Data.csv", "stats.csv"]
    orbit_num = []
    frame     = []
    chi_2     = []
    L_2       = []
    L_inf     = []


    for i in range(len(RUN_NAME)):
        for j in range(len(ANALYSIS_NAME)):
            file_path = "Plots/" + RUN_NAME[i] + "/" + ANALYSIS_NAME[j] + "/" + file_name[1]

            columns = read_csv(file_path)

            orbit_num.append(columns['Orbit Number'])
            #frame.append(columns['Frame'])
            chi_2.append(columns['Chi Squared'])
            L_2.append(columns['L2'])
            L_inf.append(columns['L Infinity'])

    return orbit_num, frame, chi_2, L_2, L_inf

def get_ecc_ang_data(RUN_NAME, ANALYSIS_NAME):
    file_name = "ecc_ang.csv"
    orbit_num = []
    frame     = []
    disk_misalign = []
    ecc_inner_avg = []
    ecc_outer_avg = []
    ang_inner_avg = []
    ang_outer_avg = []
    

    for i in range(len(RUN_NAME)):
        for j in range(len(ANALYSIS_NAME)):
            file_path = "Plots/" + RUN_NAME[i] + "/" + ANALYSIS_NAME[j] + "/" + file_name

            columns = read_csv(file_path)

            orbit_num.append(columns['Orbit Number'])
            frame.append(columns['Frame'])
            disk_misalign.append(columns['Disk Misalignment'])
            ecc_inner_avg.append(columns['Ecc Inner'])
            ecc_outer_avg.append(columns['Ecc Outer'])
            ang_inner_avg.append(columns['Angle Inner'])
            ang_outer_avg.append(columns['Angle Outer'])

    return orbit_num, frame, disk_misalign, ecc_outer_avg, ecc_inner_avg, ang_inner_avg, ang_outer_avg


def load_spec_data():
    #loads the spectrum data and thermal data
    file = open("Data/high_J.txt","r+")
    vel_j  = []
    ein_j  = []
    eng_j  = []
    g1j    = []
    g0j    = []

    for line in file:
        co_data = line.split()
        vel_j.append(float(co_data[1]))
        ein_j.append(float(co_data[2]))
        eng_j.append(float(co_data[3]))
        g1j.append(float(co_data[8]))
        g0j.append(float(co_data[9]))
    file.close()
    print("Spectrum Data Loaded")

    file = open("Data/Partfun_12C16O.txt","r+")
    temp  = []
    part  = []

    for line in file:
        part_data = line.split()
        temp.append(float(part_data[0]))
        part.append(float(part_data[1]))
    file.close()
    return vel_j, ein_j, eng_j, g1j, g0j, temp, part


def process_files(LABEL, info):
    PATH, DIM, IMAX, JMAX, KMAX = info
    time_spacing = []
    time_spacing_arrays = []

    last_processed_time = None
    new_array_flag = False

    for k in range(0,99999):
        if DIM=='2D':
            exist = r.check_file_exist_2D(PATH, IMAX, JMAX, LABEL, k)
            full_path = PATH+"binary_"+str(IMAX)+"x"+str(JMAX)+"_"+LABEL+"_"+r.frame_num(k)

        if DIM=='3D':
            exist = r.check_file_exist_3D(PATH, IMAX, JMAX, KMAX, LABEL, k)
            full_path = PATH+"binary_"+str(IMAX)+"x"+str(JMAX)+"x"+str(KMAX)+"_"+LABEL+"_"+r.frame_num(k)

        if exist:
            current_time = os.path.getctime(full_path)
            creation_date_time = datetime.fromtimestamp(current_time)
            
            if last_processed_time:
                time_difference = current_time - last_processed_time
                if time_difference > 24 * 3600:
                    time_spacing_arrays.append(time_spacing)
                    time_spacing = []
                    new_array_flag = True
                else:
                    time_spacing.append(time_difference)
            
            last_processed_time = current_time
            
            if new_array_flag:
                time_spacing.append(current_time)
                new_array_flag = False
        else:
            print("Orbit frame does not exists: " +str(k))
            orbits = np.linspace(0,k,k+1)
            break
            
    time_spacing_arrays.append(time_spacing)
    return time_spacing_arrays, orbits


def load_circular_2d():
    #circular slab model
    x_len   = 1296
    y_len   = 2160
    lengths = x_len, y_len

    radius = np.zeros((x_len,y_len))
    theta  = np.zeros((x_len,y_len))
    
    rad = np.linspace(RADIUS_IN, RADIUS_OUT, x_len)
    the = np.linspace(0.0, 2*np.pi, y_len)

    for i in range(0,x_len):
        for j in range(0,y_len):
            radius[:,j] = rad
            theta[i,:]  = the
    radius = np.transpose(radius)
    theta  = np.transpose(theta)
    radius = r.cell_center_2D(radius)
    theta  = r.cell_center_2D(theta)
    grids  = radius, theta
    
    temp    = T0 * (radius/RADIUS_IN)**ALPHA
    density = N0 * (radius/RADIUS_IN)**BETA
    vel     = np.sqrt(G * MASS_STAR / (radius*AU_CM)) * np.sin(theta) * np.sin(INCLIN) / 100000
    data    = density, temp, vel

    return grids, lengths, data


def load_PEnGUIn_2d(LABEL, IMAX, JMAX, k):
    #loads 2d PEnGUIN data frame
    exist = r.check_file_exist_2D(LOAD_PATH, IMAX, JMAX, LABEL, k)

    if exist:
        data = r.load_2D_data(LOAD_PATH, IMAX, JMAX, LABEL, k)
        #time, x, y, density, pressure, x velocity, y velocity
        #0,    1, 2, 3,       4,        5,          6,

        time = np.round(data[0], 2)
        print("Time: " + str(time) + " and Orbit # " +str(k * orbit_interval))

        x_len = len(data[1])
        y_len = len(data[2])
        lengths = x_len, y_len

        radius = np.zeros((x_len,y_len))
        theta  = np.zeros((x_len,y_len))
            
        for i in range(0,x_len):
            for j in range(0,y_len):
                radius[i,j] = data[1][i]
                theta[i,j]  = data[2][j]
        radius = np.transpose(radius)
        theta = np.transpose(theta)
        radius = r.cell_center_2D(radius)
        theta  = r.cell_center_2D(theta)

        grids = radius, theta

        density  = data[3]
        pressure = data[4]
        r_vel    = data[5]
        t_vel    = data[6]

        r_vel, t_vel  = denorm_vel_2D(r_vel, t_vel)

        temp = np.round(((density / pressure)) * (G * MASS_STAR / AU_CM) * 2.3 / R )
        density = density
        data = time, pressure, density, temp, r_vel, t_vel
        return grids, lengths, data
    else:
        grids   = 0, 0
        lengths = 0, 0
        data    = 0, 0, 0, 0, 0, 0
        return grids, lengths, data

                        
def load_PEnGUIn_3d(LABEL, IMAX, JMAX, KMAX, k):
    #loads 3d PEnGUIN data frame
    exist = r.check_file_exist_3D(LOAD_PATH, IMAX, JMAX, KMAX, LABEL, k)

    if exist:
        data = r.load_3D_data(LOAD_PATH, IMAX, JMAX, KMAX, LABEL, k)
        #time, x, y, z, density, pressure, x velocity, y velocity, z velocity
        #0,    1, 2, 3, 4,       5,        6,          7,          8

        time = np.round(data[0], 2)
        print("Time: " + str(time) + " and Orbit # " +str(k * orbit_interval))

        x_len = len(data[1])
        y_len = len(data[2])
        z_len = len(data[3])
        lengths = x_len, y_len, z_len

        radius = np.zeros((x_len,y_len,z_len))
        theta  = np.zeros((x_len,y_len,z_len))
        phi    = np.zeros((x_len,y_len,z_len))

        for i in range(0,x_len):
            for j in range(0,y_len):
                for k in range(0,z_len):
                    radius[i,j,k] = data[1][i]
                    theta[i,j,k]  = data[2][j]
                    phi[i,j,k]    = data[3][j]

        grids = radius, theta, phi

        density  = data[4]
        pressure = data[5]
        r_vel    = data[6]
        t_vel    = data[7]
        p_vel    = data[8]

        r_vel, t_vel, p_vel  = denorm_vel_3D(r_vel, t_vel, p_vel)

        temp = np.round(((density / pressure)) * (G * MASS_STAR / AU_CM) * 2.3 / R )
        density = density
        data = time, pressure, density, temp, r_vel, t_vel, p_vel
        return grids, lengths, data


#-------------------------------------------------------------------------------------------------------
#saving data -------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def save_line_to_csv(data, orbit_num, frame, file_path):
    column_names = np.array(["Velocity Col " + str(i) + " [km/s]" for i in range(len(data))])

    # Check if the file already exists
    file_exists = os.path.isfile(file_path)

    with open(file_path, 'a', newline='') as csvfile:
        fieldnames = ['Orbit Number', 'Frame'] + column_names.tolist()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # If the file doesn't exist, write the header
        if not file_exists:
            writer.writeheader()

        data_dict = {'Orbit Number': orbit_num, 'Frame': frame, **{name: value for name, value in zip(column_names, data)}}
        writer.writerow(data_dict)
    return


def save_stats_to_csv(constants, file_path):
    # Check if the file already exists
    file_exists = os.path.isfile(file_path)

    # If the file doesn't exist, create it and write the header
    with open(file_path, 'a', newline='') as csvfile:
        fieldnames = ['Orbit Number', 'Frame', 'STD', 'Variance', 'Chi Squared', 'L2', 'L Infinity']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # If the file doesn't exist, write the header
        if not file_exists:
            writer.writeheader()

        constants_dict = {fieldnames[i]: constants[i] for i in range(min(len(fieldnames), len(constants)))}
        writer.writerow(constants_dict)
        return

    
def save_ecc_ang_to_csv(constants, file_path):
    # Check if the file already exists
    file_exists = os.path.isfile(file_path)

    with open(file_path, 'a', newline='') as csvfile:
        fieldnames = ['Orbit Number', 'Frame', 'Disk Misalignment', 'Ecc Inner', 'Ecc Outer', 'Angle Inner', 'Angle Outer']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # If the file doesn't exist, write the header
        if not file_exists:
            writer.writeheader()

        constants_dict = {fieldnames[i]: constants[i] for i in range(min(len(fieldnames), len(constants)))}
        writer.writerow(constants_dict)
        return
    
#-------------------------------------------------------------------------------------------------------
#file management ---------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def dir_check(run_name, analysis_name, file_name):
    #checks if directory exists
    if os.path.isdir("Plots/" + str(run_name) + "/" + str(analysis_name) + "/"):
        print("Dir exists!")
    #makes directory if not
    else:
        if os.path.isdir("Plots/" + str(run_name) + "/"):
            print("Dir analysis exists!")
        else:
            os.mkdir("Plots/" + str(run_name) + "/")
        os.mkdir("Plots/" + str(run_name) + "/" + str(analysis_name) + "/")

        for i in range(len(file_name)):
            #makes all the sub directories
            os.mkdir("Plots/" + str(run_name) + "/" + str(analysis_name) + "/" + str(file_name[i]) + "/")
        print("Dir made!")
    return    


def check_file(PATH, IMAX, JMAX, LABEL, FRAME_START, FRAME_END, DIM):
    #checks if file exists
    for k in range(FRAME_START, FRAME_END):
        if DIM=='2D':
            exist = r.check_file_exist_2D(PATH, IMAX, JMAX, LABEL, k)

        if DIM=='3D':
            exist = r.check_file_exist_3D(PATH, IMAX, JMAX, KMAX, LABEL, k)

        if exist:
            continue
        else:
            print("Orbit frame does not exists: " +str(k * orbit_interval))
            return k
        
        
def check_plots(PATH, DIR_NAME, FILE_NAME, FRAME_START, FRAME_END):
    #checks if plot exists
    for k in range(FRAME_START, FRAME_END):
        exist = os.path.isfile("Plots/" + PATH + "/" + DIR_NAME + "/" + FILE_NAME + str(k * orbit_interval) + ".png")

        if exist:
            continue
        else:
            print("Plot frame does not exists: " +str(k * orbit_interval))
            return k


#-------------------------------------------------------------------------------------------------------
#velocity calculations ---------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def denorm_vel_2D(r_vel, t_vel):
    #denormalizes 2d velocity data to km/s
    r_vel = r_vel * np.sqrt(G * MASS_STAR / (AU_CM))/100000
    t_vel = t_vel * np.sqrt(G * MASS_STAR / (AU_CM))/100000
    return r_vel, t_vel


def denorm_vel_3D(r_vel, t_vel, p_vel):
    #denormalizes 3d velocity data to km/s
    r_vel = r_vel * np.sqrt(G * MASS_STAR / (AU_CM))/100000
    t_vel = t_vel * np.sqrt(G * MASS_STAR / (AU_CM))/100000
    p_vel = p_vel * np.sqrt(G * MASS_STAR / (AU_CM))/100000
    return r_vel, t_vel, p_vel


def incvel_2D(r_vel, t_vel, theta, run_type, roll_num):
    #inclines 2d velocity data to viewing geometry
    if run_type=="run_rotate":
        r_vel = np.roll(r_vel, roll_num, axis=0)
        t_vel = np.roll(t_vel, roll_num, axis=0)

    vel_y = (r_vel  * np.sin(theta) + (t_vel  * np.cos(theta)))
    vzp   =  vel_y * np.sin(INCLIN)
    return vzp


def incvel_3D(r_vel, t_vel, p_vel, theta, run_type, roll_num):
    #inclines 2d velocity data to viewing geometry
    if run_type=="run_rotate":
        r_vel = np.roll(r_vel, roll_num, axis=0)
        t_vel = np.roll(t_vel, roll_num, axis=0)
        p_vel = np.roll(p_vel, roll_num, axis=0)

    vel_y = (r_vel * np.sin(theta) * np.sin(np.pi/2) +
            (p_vel * np.sin(theta) * np.cos(np.pi/2) + t_vel * np.cos(theta) * np.sin(np.pi/2)))
    vel_z = (r_vel * np.cos(np.pi/2) - p_vel * np.sin(np.pi/2))
    vzp   =  vel_y * np.sin(INCLIN) + vel_z * np.cos(INCLIN)
    return r_vel, t_vel, p_vel, vzp


#-------------------------------------------------------------------------------------------------------
#flux and line density calculations --------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def flux_cal(density, vel_j, ein_j, eng_j, g1j, g0j, temp, part_func):
    #calculates the flux of the disk for given transition lines
    constant1 = 1 / (8 * B_CGS * (np.pi)**(3./2.))
    x, y      = np.shape(density)
    optical_t = np.zeros((len(g0j), x, y))
    flux_dens = np.zeros((x, y))

    for i in range(len(g0j)):
        #loops through the different transitions
        opt_depth        = (density / part_func)  * np.exp(-eng_j[i] / (BOLTZ * temp)) * (ein_j[i] * g0j[i]**2 / (g1j[i] * vel_j[i]**3)) * constant1
        optical_t[i,:,:] = opt_depth.astype(float)
        planck_func      = (2 * PLANCK * C * vel_j[i]**3) / (np.exp((PLANCK * C * vel_j[i])/(BOLTZ_CGS * temp)) - 1)
        flux_dens       += (1-np.exp(-opt_depth.astype(float))) * vel_j[i]/ C * planck_func
    return flux_dens


def flux_cal_norm(density, vel_j, ein_j, eng_j, g1j, g0j, temp, part_func):
    #calculates the flux of the disk for given normalized transition lines
    constant1 = 1 / (8 * B_CGS * (np.pi)**(3./2.))
    x, y      = np.shape(density)
    optical_t = np.zeros((len(g0j), x, y))
    flux_dens = np.zeros((np.shape(density)))

    for i in range(len(g0j)):
        #loops through the different transitions
        opt_depth        = (density / part_func)  * np.exp(-eng_j[i] / (BOLTZ * temp)) * (ein_j[i] * g0j[i]**2 / (g1j[i] * vel_j[i]**3)) * constant1
        optical_t[i,:,:] = opt_depth.astype(float)
        planck_func      = (2 * PLANCK * C * vel_j[i]**3) / (np.exp((PLANCK * C * vel_j[i])/(BOLTZ_CGS * temp)) - 1)
        line_grid        = (1-np.exp(-opt_depth.astype(float))) * vel_j[i]/ C * planck_func
        flux_dens       += line_grid / np.max(line_grid)
    return flux_dens, optical_t


def int_j_cal(flux_dens, radius, k_vel, vzp):
    #line density calculation
    int_j = np.zeros(np.shape(k_vel))
    constant2 = 1 / (B_CGS * np.sqrt(np.pi))

    flux_radius = flux_dens * radius**2 * AU_KM
    for i in range(len(k_vel)):
        #loops through the test speeds
        int_j[i] = np.sum(flux_radius * constant2 * np.exp(-(k_vel[i] - vzp)**2 / B**2))
    return int_j


def part_int_j_cal(flux_dens, radius, k_vel, vzp):
    vel  = vzp[np.logical_not(np.isnan(vzp))]
    flux = flux_dens[np.logical_not(np.isnan(vzp))]
    rad  = radius[np.logical_not(np.isnan(vzp))]

    #line density calculation
    int_j = np.zeros(np.shape(k_vel))
    constant2 = 1 / (B_CGS * np.sqrt(np.pi))

    flux_radius = flux * rad**2 * AU_KM
    for i in range(len(k_vel)):
        #loops through the test speeds
        int_j[i] = np.sum(flux_radius * constant2 * np.exp(-(k_vel[i] - vel)**2 / B**2))
        #think about flux conservation
    return int_j

def vzp_masking(min_value, max_value, array):
    #creates mask
    mask = (array < min_value) | (array > max_value)

    array[mask] = None
    
    return array


#-------------------------------------------------------------------------------------------------------
#thermal calculations ----------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def part_func_1D(temp, part, x_len):
    #generates 1d partion function grid from 1d radial temperatures
    part_func = np.zeros((np.shape(temp)))
    for i in range(0, x_len-1):
        if temp[0,i] > 8999.0:
            part_func[:,i] = part[-1]
        else:
            part_func[:,i] = part[int(temp[0,i])]
    return part_func


def part_func_2D(temp, part, x_len, y_len):
    #generates 2d partion function grid from 1d radial, azimuthal temperatures
    part_func = np.zeros((np.shape(temp)))
    print(np.shape(temp))
    for i in range(0,y_len-1):
        for j in range(0,len(temp[0,:])-1):
            #for j in range(0,x_len-1):
            if temp[i,j] > 8999.0:
                part_func[i,j] = part[-1]
            else:
                part_func[i,j] = part[int(temp[i,j])]
    return part_func


#-------------------------------------------------------------------------------------------------------
#data processing ---------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def remove_spike(int_j):
    #removes atmospheric spike in observation, synthetic, and velocity
    spots_to_remove = [296,297,298,299,300,301,302,303,304]
    
    new_j = int_j
    new_j = np.delete(new_j, spots_to_remove)
    return new_j


def equal_area_data(vzp, new_j, new_obs, new_k):
    #trims and equal area normalizes the data
    min_x, max_x = 280-int(abs(np.floor(np.min(vzp)))), 311+int(np.ceil(np.max(vzp)))

    norm_obs = new_obs[min_x:max_x]
    norm_j   = new_j[min_x:max_x]
    norm_k   = new_k[min_x:max_x]
    
    #equal area normalization
    area1 = trap_area(norm_k, norm_obs)
    area2 = trap_area(norm_k, norm_j)
    
    norm_obs = norm_obs / area1
    norm_j   = norm_j / area2
    return norm_obs, norm_j, norm_k, area1, area2


#-------------------------------------------------------------------------------------------------------
#stat functions ----------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def correct_orbit_num(orbit_num, correction_num):
    #fixes orbit num with rotational analysis
    corrected_orbit_num = []
    for i in range(len(correction_num)):
        corrected_orbit_num.append((orbit_num[i]-correction_num[i])*2*np.pi/180)
    return corrected_orbit_num

def correct_frame(frame, correction_num):
    #fixes frame num with rotational analysis
    corrected_frame = []
    for i in range(len(correction_num)):
        corrected_frame.append((frame[i]-correction_num[i])*2*np.pi/180)
    return corrected_frame


#-------------------------------------------------------------------------------------------------------
#scaling functions -------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def scale_data(radius, r_vel, t_vel):
    #scales the data from code units to physical units
    radius = radius * PLANET_SCALING
    r_vel  = r_vel * ORBIT_FREQ * PLANET_SCALING
    t_vel  = t_vel * ORBIT_FREQ * PLANET_SCALING
    return radius, r_vel, t_vel


#-------------------------------------------------------------------------------------------------------
#misc functions ----------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def cal_mass_disk(density):
    #returns mass estimate of disk in Mj
    dens_mass = np.sum(np.trapz(density)) * MOL * H_MASS * H_CO_RATIO / MASS_JUP
    return dens_mass


def break_data(breaks, data, params):
    #breaks data in segments
    p_dens, p_temp, r_vel, t_vel, radius, theta, x_len, y_len = data
    RUN_TYPE, roll_num = params
    
    vel_j, ein_j, eng_j, g1j, g0j, tempature, part = load_spec_data()
    
    x, y  = np.shape(radius)
    int_j = np.zeros((np.shape(k_vel)))

    # Initialize lists to store the segmented arrays
    b_density = []
    b_temp = []
    b_radius = []
    b_theta = []
    b_r_vel = []
    b_t_vel = []
    b_x_len = []

    b_vzp = []
    b_part_func = []
    b_flux_dens = []
    b_optical_t = []
    b_int_j = []

    for i in range(len(breaks)-1):
        start_x = breaks[i]
        end_x = breaks[i+1]

        subarray_shape = (x, end_x - start_x)

        # Create arrays for each segment
        b_density.append(p_dens[:, start_x:end_x])
        b_temp.append(p_temp[:, start_x:end_x])
        b_r_vel.append(r_vel[:, start_x:end_x])
        b_t_vel.append(t_vel[:, start_x:end_x])
        b_radius.append(radius[:, start_x:end_x])
        b_theta.append(theta[:, start_x:end_x])
        b_x_len.append(int(end_x - start_x + 1))

        # Calculate other values for the segment
        b_vzp.append(incvel_2D(b_r_vel[i], b_t_vel[i], b_theta[i], RUN_TYPE, roll_num))
        b_part_func.append(part_func_2D(b_temp[i], part, int(b_x_len[i]), int(y_len)))
        b_flux_dens.append(flux_cal(b_density[i], vel_j, ein_j, eng_j, g1j, g0j, b_temp[i], b_part_func[i]))
        b_int_j.append(int_j_cal(b_flux_dens[i], b_radius[i], k_vel, b_vzp[i]))
        #print("Break " + str(i) + "/" + str(len(breaks)-1) + " finished!")
        
        int_j += b_int_j[i]
    return int_j, b_int_j


#-------------------------------------------------------------------------------------------------------
#ecc calculations --------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
def trap_area(x, curve):
    #returns area of 1d curve via trapozoidal method
    y0 = curve
    y1 = np.roll(curve, 1, axis=0)
    
    x0 = x[0]
    x1 = x[1]
    
    values = 0.5 * (y0 + y1) * (x1 - x0)
    
    area = np.sum(values, axis=0)
    return area

def ecc_calculation(theta, r_vel, t_vel):
    #returns the fourier analysis ecc and perapsis
    curve = np.cos(theta) * r_vel 
    area  = trap_area(theta, curve)
    A     = ECC_NORM * area
    
    curve = np.sin(theta) * r_vel 
    area  = trap_area(theta, curve)
    B     = ECC_NORM * area

    val   = np.sqrt(A**2 + B**2)
    t_avg = np.average(t_vel, axis =0)
    ecc   = val / t_avg
    ang   = np.arctan2(B,A)
    
    #pi/2 angle correction for sine, cosine 
    ang = ang - np.pi/2

    return ecc, ang

def ecc_ang_avgs(radius, ecc, ang, density):
    #returns the density weighted average ecc and perapsis
    inner_break = np.abs(radius[:] - INNER).argmin() #inner edge to 80% of planet orbit
    outer_break = np.abs(radius[:] - OUTER).argmin() #120% of planet orbit to outer edge
    
    density = np.average(density, axis=0)

    inner_dr = radius[1:inner_break+1] - radius[:inner_break]
    outer_dr = radius[outer_break:] - radius[outer_break-1:-1]

    ecc_inner_avg = np.round(np.sum(ecc[:inner_break]*inner_dr*density[:inner_break])/np.sum(inner_dr*density[:inner_break]),4)
    ecc_outer_avg = np.round(np.sum(ecc[outer_break:]*outer_dr*density[outer_break:])/np.sum(outer_dr*density[outer_break:]),4)
    ecc_avg       = [ecc_inner_avg, ecc_outer_avg]

    #need to check if 0/2pi looping values is still correct in the average
    ang_inner_avg = np.round(np.sum(ang[:inner_break]*inner_dr*density[:inner_break])/np.sum(inner_dr*density[:inner_break]),4) % (2 * np.pi)
    ang_outer_avg = np.round(np.sum(ang[outer_break:]*outer_dr*density[outer_break:])/np.sum(outer_dr*density[outer_break:]),4) % (2 * np.pi)
    ang_avg       = [ang_inner_avg, ang_outer_avg]

    disk_misalign = np.round((ang_inner_avg - ang_outer_avg) % (2 * np.pi), 4)
    ang = ang % (2 * np.pi)
    
    return disk_misalign, ecc_avg, ang_avg, ang


def moving_slope(data, window_size):
    #need to check if 0/2pi looping values is still correct in the average
    slopes = []
    for i in range(len(data) - window_size + 1):
        x = np.arange(i, i + window_size)
        slope, _, _, _, _ = linregress(x, data[i:i+window_size])
        slopes.append(slope)
    return np.array(slopes)

def circular_statistics(data, window_size=3):
    #need to check if 0/2pi looping values is still correct in the average
    variance      = []
    std_deviation = []
    mean_angle    = []
    slope         = []
    
    for i in range(3):
        variance.append(np.round(circvar(data[i]),4))
        std_deviation.append(np.round(circstd(data[i]),4))
        mean_angle.append(np.round(circmean(data[i]),4))
        slope.append(np.round(moving_slope(data[i], window_size),4))

    return variance, std_deviation, mean_angle, slope