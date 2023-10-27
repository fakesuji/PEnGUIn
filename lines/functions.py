"""
Author:  Cory Padgett
Advisor: Dr. Jeffrey Fung
Email:   cpadge4@clemson.edu
"""

#imports
import os
import numpy as np
import read_penguin as r


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


#system constants
MASS_STAR   = 1.02 * MASS_SOL   #g
MASS_PLANET = 11.6 * MASS_JUP   #g
INCLIN  = np.radians(71.0)      #radians
T0      = 2770           #K
ALPHA   = -0.33          #
N0      = 1.0E19
BETA    = -2.0
R_IN    = 0.048          #au
R_OUT   = 0.898          #au
B       = 8.8            #km/s
B_CGS   = B * 100000     #cm/s
INNER   = 0.13
OUTER   = 0.19
DENS_FACTOR = 9E19


#arrays
k_vel      = np.linspace(-300,300,601)


def load_circular_2d():
    x_len   = 1296
    y_len   = 2160
    lengths = x_len, y_len

    radius = np.zeros((x_len,y_len))
    theta  = np.zeros((x_len,y_len))
    
    rad = np.linspace(R_IN, R_OUT, x_len)
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
    
    temp    = T0 * (radius/R_IN)**ALPHA
    density = N0 * (radius/R_IN)**BETA
    vel     = np.sqrt(G * MASS_STAR / (radius*AU_CM)) * np.sin(theta) * np.sin(INCLIN) / 100000
    data    = density, temp, vel

    return grids, lengths, data


def load_PEnGUIn_2d(LABEL, k):
    LOAD_PATH  = "/scratch/cpadge4/"
    SAVE_PATH  = "/scratch1/cpadge4/"

    IMAX  = 1200
    JMAX  = 2160

    exist = r.check_file_exist_2D(LOAD_PATH, IMAX, JMAX, LABEL, k)

    if exist:
        data = r.load_2D_data(LOAD_PATH, IMAX, JMAX, LABEL, k)
        #time, x, y, density, pressure, x velocity, y velocity
        #0,    1, 2, 3,       4,        5,          6,

        time = np.round(data[0], 2)
        print("Time: " + str(time) + " and Orbit # " +str(k))

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

        pressure = data[3]
        density  = data[4]
        r_vel    = data[5]
        t_vel    = data[6]

        r_vel, t_vel  = denorm_vel_2D(r_vel, t_vel)

        temp = np.round(((density / pressure)) * (G * MASS_STAR / AU_CM) * 2.3 / R )
        density = density  * DENS_FACTOR
        data = time, pressure, density, temp, r_vel, t_vel
        return grids, lengths, data

                        
def load_PEnGUIn_3d(LABEL, k, DENS_FACTOR):
    LOAD_PATH  = "/scratch/cpadge4/"
    SAVE_PATH  = "/scratch1/cpadge4/"

    IMAX  = 624
    JMAX  = 1152
    KMAX  = 24

    exist = r.check_file_exist_3D(LOAD_PATH, IMAX, JMAX, KMAX, LABEL, k)

    if exist:
        data = r.load_3D_data(LOAD_PATH, IMAX, JMAX, KMAX, LABEL, k)
        #time, x, y, z, density, pressure, x velocity, y velocity, z velocity
        #0,    1, 2, 3, 4,       5,        6,          7,          8

        time = np.round(data[0], 2)
        print("Time: " + str(time) + " and Orbit # " +str(k))

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

        pressure = data[4]
        density  = data[5]
        r_vel    = data[6]
        t_vel    = data[7]
        p_vel    = data[8]

        r_vel, t_vel, p_vel  = denorm_vel_3D(r_vel, t_vel, p_vel)

        temp = np.round(((density / pressure)) * (G * MASS_STAR / AU_CM) * 2.3 / R )
        density = density  * DENS_FACTOR
        data = time, pressure, density, temp, r_vel, t_vel, p_vel
        return grids, lengths, data


def cal_mass_disk(density):
    dens_mass = np.sum(np.trapz(density)) * MOL * H_MASS * H_CO_RATIO / MASS_JUP
    return dens_mass
    

def denorm_vel_2D(r_vel, t_vel):
    r_vel = r_vel * np.sqrt(G * MASS_STAR / (AU_CM))/100000
    t_vel = t_vel * np.sqrt(G * MASS_STAR / (AU_CM))/100000
    return r_vel, t_vel


def denorm_vel_3D(r_vel, t_vel, p_vel):
    r_vel = r_vel * np.sqrt(G * MASS_STAR / (AU_CM))/100000
    t_vel = t_vel * np.sqrt(G * MASS_STAR / (AU_CM))/100000
    p_vel = p_vel * np.sqrt(G * MASS_STAR / (AU_CM))/100000
    return r_vel, t_vel, p_vel


def incvel_2D(r_vel, t_vel, theta, run_type, roll_num):
    if run_type=="run_rotate":
        r_vel = np.roll(r_vel, roll_num, axis=1)
        t_vel = np.roll(t_vel, roll_num, axis=1)

    vel_y  = (r_vel * np.sin(theta) + (t_vel * np.cos(theta)))
    vzp = vel_y * np.sin(INCLIN)
    return vzp


def incvel_3D(r_vel, t_vel, p_vel, theta, run_type, roll_num):
    if run_type=="run_rotate":
        r_vel = np.roll(r_vel, roll_num, axis=1)
        t_vel = np.roll(t_vel, roll_num, axis=1)
        p_vel = np.roll(p_vel, roll_num, axis=1)

    vel_y  = (r_vel * np.sin(theta) * np.sin(np.pi/2) +
           (p_vel * np.sin(theta) * np.cos(np.pi/2) + t_vel * np.cos(theta) * np.sin(np.pi/2)))
    vel_z  = (r_vel * np.cos(np.pi/2) - p_vel * np.sin(np.pi/2))
    vzp = vel_y * np.sin(INCLIN) + vel_z * np.cos(INCLIN)
    return r_vel, t_vel, p_vel, vzp


def flux_cal(density, vel_j, ein_j, eng_j, g1j, g0j, temp, part_func):
    constant1 = 1 / (8 * B_CGS * (np.pi)**(3./2.))
    x,y = np.shape(density)
    optical_t = np.zeros((len(g0j),x,y))
    flux_dens = np.zeros((np.shape(density)))

    for i in range(len(g0j)):
        opt_depth        = (density / part_func)  * np.exp(-eng_j[i] / (BOLTZ * temp)) * (ein_j[i] * g0j[i]**2 / (g1j[i] * vel_j[i]**3)) * constant1
        optical_t[i,:,:] = opt_depth.astype(float)
        planck_func      = (2 * PLANCK * C * vel_j[i]**3) / (np.exp((PLANCK * C * vel_j[i])/(BOLTZ_CGS * temp)) - 1)
        flux_dens       += (1-np.exp(-opt_depth.astype(float))) * vel_j[i]/ C * planck_func
    return flux_dens, optical_t


def flux_cal_norm(density, vel_j, ein_j, eng_j, g1j, g0j, temp, part_func):
    constant1 = 1 / (8 * B_CGS * (np.pi)**(3./2.))
    x,y = np.shape(density)
    optical_t = np.zeros((len(g0j),x,y))
    flux_dens = np.zeros((np.shape(density)))

    for i in range(len(g0j)):
        opt_depth        = (density / part_func)  * np.exp(-eng_j[i] / (BOLTZ * temp)) * (ein_j[i] * g0j[i]**2 / (g1j[i] * vel_j[i]**3)) * constant1
        optical_t[i,:,:] = opt_depth.astype(float)
        planck_func      = (2 * PLANCK * C * vel_j[i]**3) / (np.exp((PLANCK * C * vel_j[i])/(BOLTZ_CGS * temp)) - 1)
        line_grid        = (1-np.exp(-opt_depth.astype(float))) * vel_j[i]/ C * planck_func
        flux_dens       += line_grid / np.max(line_grid)
    return flux_dens, optical_t


def int_j_cal(flux_dens, radius, k_vel, vzp):
    int_j = np.zeros(np.shape(k_vel))
    constant2 = 1 / (B_CGS * np.sqrt(np.pi))

    flux_radius = flux_dens * radius**2 * AU_KM
    for i in range(len(k_vel)):
        int_j[i] = np.sum(flux_radius * constant2 * np.exp(-(k_vel[i] - vzp)**2 / B**2))
    return int_j


def part_func_1D(temp, part, y_len):
    part_func = np.zeros((np.shape(temp)))
    for i in range(0, y_len-1):
        if temp[i,0] > 8999.0:
            part_func[i,:] = part[-1]
        else:
            part_func[i,:] = part[int(temp[i,0])]
    return part_func


def part_func_2D(temp, part, x_len, y_len):
    part_func = np.zeros((np.shape(temp)))
    for i in range(0,y_len-1):
        for j in range(0,x_len-1):
            if temp[i,j] > 8999.0:
                part_func[i,j] = part[-1]
            else:
                part_func[i,j] = part[int(temp[i,j])]
    return part_func


def dir_check(run_name, analysis_name, file_name, dim):
    #checks if directory exists
    if os.path.isdir("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/"):
        print("Dir exists!")
    #makes directory if not
    else:
        if os.path.isdir("Plots/run_" + str(run_name) + "/"):
            print("Dir analysis exists!")
        else:
            os.mkdir("Plots/run_" + str(run_name) + "/")
        os.mkdir("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/")

        if dim=='2D':
            for i in range(2,14):
                #makes all the sub directories
                os.mkdir("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + str(file_name[i]) + "/")
            print("Dir made!")
        if dim=='3D':
            for i in range(3,14):
                #makes all the sub directories
                os.mkdir("Plots/run_" + str(run_name) + "/" + str(analysis_name) + "/" + str(file_name[i]) + "/")
            print("Dir made!")
    return


def load_spec_data():
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


def make_dim_match(x, y, mapible):
    x_shape = np.shape(x)
    y_shape = np.shape(y)
    m_shape = np.shape(mapible)
    
    if m_shape == (x_shape, y_shape):
        print("Dim Match!")
    else:
        xm_len = m_shape[0]
        ym_len = m_shape[1]
        
        x = x[:xm_len]
        y = y[:ym_len]
    return x, y