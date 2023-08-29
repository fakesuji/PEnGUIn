#imports
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import read_penguin as r

#general constants
h      = 6.626068 * 10**(-27)      #erg s       Planck’s Constant
kb     = 0.6950348                 #cm^-1 K^-1  Boltzman’s Constant
kb_cgs = 1.380649 * 10**(-16)      #erg K^-1  Boltzman’s Constant
c      = 2.997925 * 10**10         #cm s^-1     Speed of light in vacuum
G      = 6.674 * 10**(-8)          #cm^3 g^−1 s^−2
R      = 8.31446261815324 * 10**7  #erg K^-1 mol^-1
mol    = 6.022 * 10*23
au_cm    = 1.496 * 10**13     #cm
au_km    = 1.496 * 10**8      #km
Mass_Sol = 1.989 * 10**33     #g
Mass_Jup = 1.898 * 10**30     #g


#system constants
Mass_Star   = 1.02 * Mass_Sol   #g
Mass_Planet = 11.6 * Mass_Jup   #g
inclin  = np.radians(71.0)      #radians
b       = 8.8            #km/s

#dir creation, indexing frames, and setting run info
run_num     = 3.1
run_type    = "run_rotate"  #"run_stationary" or "run_rotate"
re_run      = "NO"          #"YES" or "NO"
orbit_num   = 1005          #Used only in run_rotate run type
frame_start = 1769          #initial frame to pull
frame_end   = 1897          #final frame -1 to pull
dim         = '2D'          #"2D" or "3D"
dens_factor = 1.0E19

name_3D = ['Time', 'X position', 'Y position', 'Z Position', 
           'Density', 'Pressure', 'X Velocity' ,'Y Velocity', 'Z Velocity', 
           'Multiple Plot', 'Projected Velocity','Ploar Density', 'Line Density']
file_name_3D = ['Time', 'X_position', 'Y_position', 'Z Position', 
                'Density', 'Pressure', 'X_Velocity' ,'Y_Velocity', 'Z_Velocity', 
                'Multi_Plot', 'VZP', 'polar_dens', 'Line_Density']

name_2D = ['Time', 'X position', 'Y position', 'Density', 
           'Pressure', 'X Velocity' ,'Y Velocity', 'Multiple Plot', 
           'Projected Velocity','Ploar Density', 'Line Density']
file_name_2D = ['Time', 'X_position', 'Y_position', 'Density', 
                'Pressure', 'X_Velocity' ,'Y_Velocity', 
                'Multi_Plot', 'VZP', 'polar_dens', 'Line_Density']

#data to be saved in text files
save_file_name = ["part_func", "flux_dens", "opt_depth", "Ij", "Vzp"]

#checks if directory exists
if os.path.isdir("Plots/run_" + str(run_num) + "/"):
    print("Dir exists!")
#makes directory if not
else:
    os.mkdir("Plots/run_" + str(run_num) + "/")
    os.mkdir("Data/run_" + str(run_num) + "/")
    for i in range(len(save_file_name)):
        os.mkdir("Data/run_" + str(run_num) + "/" + save_file_name[i] + "/")
        
    if dim=='2D':
        for i in range(3,11):
            #makes all the sub directories
            os.mkdir("Plots/run_" + str(run_num) + "/" + str(file_name_2D[i]) + "/")
        print("Dir made!")
    if dim=='3D':
        for i in range(4,13):
            #makes all the sub directories
            os.mkdir("Plots/run_" + str(run_num) + "/" + str(file_name_3D[i]) + "/")
        print("Dir made!")
    
    
#----------------------------------------------------------------------------------
#Spectrum Data
#----------------------------------------------------------------------------------
file = open("Data/12CO_data.txt","r+")
vj  = []
Aj  = []
Ej  = []
g1j = []
gj  = []

for line in file:
    CO_data = line.split()
    vj.append(float(CO_data[1]))
    Aj.append(float(CO_data[2]))
    Ej.append(float(CO_data[3]))
    g1j.append(float(CO_data[8]))
    gj.append(float(CO_data[9]))
file.close()


file = open("Data/Partfun_12C18O.txt","r+")
temp  = []
part  = []

for line in file:
    part_data = line.split()
    temp.append(float(part_data[0]))
    part.append(float(part_data[1]))
file.close()
print("Spectrum Data Loaded")


#--------------------------------------------------------------------------------
#load data
#--------------------------------------------------------------------------------
path  = "/scratch/cpadge4/"
if dim=='3D':
    imax  = 624
    jmax  = 1152
    kmax  = 24
    
if dim=='2D':
    imax  = 720
    jmax  = 1536
    
label = "h33_1p3J_a-10_OA_PPM4"

#binary_720x1128x24_h35_1p2J_a-20_b-20_OA_PPM4_01234 - example data file name

for k in range(frame_start,frame_end):
    if dim=='2D':
        exist = r.check_file_exist_2D(path, imax, jmax, label, k)
    
    if dim=='3D':
        exist = r.check_file_exist_3D(path, imax, jmax, kmax, label, k)
        
    if exist:
        print(k)
        if dim=='2D':
            data = r.load_2D_data(path, imax, jmax, label, k)
            time = np.round(data[0], 3)
            print(time)
            for i in range(3,7):
                fig, ax = plt.subplots()
                if i <= 4:
                    plt.pcolormesh(data[1], data[2], np.log10(data[i]))
                else:
                    plt.pcolormesh(data[1], data[2], data[i])

                plt.colorbar()
                plt.xlabel("Radius [au]")
                plt.ylabel("Theta [rads]")
                plt.title(str(name_2D[i]) + " Plot at Time " + str(time) + " and Orbit " + str(k))
                plt.savefig("Plots/run_" + str(run_num) + "/" + str(file_name_2D[i]) + "/" + str(file_name_2D[i]) + '_' + str(k) + ".png")
                plt.close()
            
        if dim=='3D':
            data = r.load_3D_data(path, imax, jmax, kmax, label, k)
            #time, x, y, z, density, pressure, x velocity, y velocity, z velocity
            time = np.round(data[0], 3)
            for i in range(4,9):
                fig, ax = plt.subplots()
                if i <= 5:
                    plt.pcolormesh(data[1], data[2], np.log10(data[i][-1,:,:]))
                else:
                    plt.pcolormesh(data[1], data[2], data[i][-1,:,:])

                plt.colorbar()
                plt.title(str(name_3D[i]) + " Plot at Time " + str(time) + ", Orbit " + str(k) + ", and Z plane at " + str(-1))
                plt.savefig("Plots/run_" + str(run_num) + "/" + str(file_name_3D[i]) + "/" + str(file_name_3D[i]) + '_' + str(k) + ".png")
                plt.close()
        
        
        #----------------------------------------------------------------------------------------------
        #initialize radius and theta and dim lengths
        #----------------------------------------------------------------------------------------------
        if dim=='2D':
            x = len(data[6][0,:])
            y = len(data[6][:,0])

        if dim=='3D':
            x = len(data[6][0,0,:])
            y = len(data[6][0,:,0])
            z = len(data[6][:,0,0])

        radius = np.zeros((x,y))
        theta  = np.zeros((x,y))


        for i in range(0,x):
            for j in range(0,y):
                radius[i,j] = data[1][i]
                theta[i,j]  = data[2][j]

        #computational constants
        r_in  = radius[0,0]  * au_cm  #cm
        r_out = radius[-1,0] * au_cm  #cm
        roll_num = 12*(k-frame_start)     #number to roll array's by for frame num
        
        inner = 0.105
        outer = 0.145

        radius_1 = np.ones((y+1,))*inner
        radius_2 = np.ones((y+1,))*outer
        
        
        #----------------------------------------------------------------------------------------------
        #make all of arrays
        #----------------------------------------------------------------------------------------------
        k_vel      = np.linspace(-300,300,601)
        part_func  = np.zeros((x, y,))
        flux_dens  = np.zeros((x, y,))
        Ij         = np.zeros((len(k_vel),))

        if dim=='2D':
            pressure = data[3][:,:]
            density = data[4][:,:]
            r_vel = np.transpose(data[5])
            t_vel = np.transpose(data[6])
                        
        if dim=='3D':
            pressure = data[4][-1,:,:]
            density = data[5][-1,:,:]
            r_vel = np.transpose(data[6][-1,:,:])
            t_vel = np.transpose(data[7][-1,:,:])
            p_vel = np.transpose(data[8][-1,:,:])
        
        
        #----------------------------------------------------------------------------------------------
        #Compute line
        #----------------------------------------------------------------------------------------------
        #----------------------------------------------------------------------------------------------
        #Flux Calculations
        #----------------------------------------------------------------------------------------------
        T = np.round(((density / pressure)) * (G * Mass_Star / au_cm) * 2.3 / R )
        density = density * dens_factor
        
        if run_type=="run_rotate":
            T       = np.roll(T, roll_num, axis=0)
            density = np.roll(density, roll_num, axis=0)

        T = np.transpose(T)
        density = np.transpose(density)
        
        if re_run == "NO":
            for i in range(0,x):
                for j in range(0,y):
                    spot = int(np.argwhere(temp==T[i,j]))
                    part_func[i,j] = part[spot]
        
        if re_run == "YES":
            part_func = np.loadtxt("Data/run_" + str(run_num) + "/data_" + str(k) + "_part_func.txt")

        constant1 = 1 / (8*b*100000*(np.pi)**(3./2.))
        constant2 = 1 / (b*100000 * np.sqrt(np.pi))

        for i in range(0,len(gj)):
            opt_depth   = (density / part_func)  * np.exp(-Ej[i] / (kb * T)) * (Aj[i] * gj[i]**2 / (g1j[i] * vj[i]**3)) * constant1
            planck_func = (2 * h * c * vj[i]**3) / (np.exp((h * c * vj[i])/(kb_cgs * T)) - 1)
            flux_dens  += (1-np.exp(-opt_depth)) * vj[i]/c * planck_func

        
        #----------------------------------------------------------------------------------------------
        #Velocity Calculations
        #----------------------------------------------------------------------------------------------
        if dim == '2D':
            r_vel = r_vel * np.sqrt(G * Mass_Star / (au_cm))/100000
            t_vel = t_vel * np.sqrt(G * Mass_Star / (au_cm))/100000
            
            if run_type=="run_rotate":
                r_vel = np.roll(r_vel, roll_num, axis=1)
                t_vel = np.roll(t_vel, roll_num, axis=1)

            Vy  = (r_vel * np.sin(theta) + (t_vel * np.cos(theta)))
            Vzp = -Vy * np.sin(inclin)
            
        if dim == '3D':
            r_vel = r_vel * np.sqrt(G * Mass_Star / (au_cm))/100000
            t_vel = t_vel * np.sqrt(G * Mass_Star / (au_cm))/100000
            p_vel = p_vel * np.sqrt(G * Mass_Star / (au_cm))/100000
            
            if run_type=="run_rotate":
                r_vel = np.roll(r_vel, roll_num, axis=1)
                t_vel = np.roll(t_vel, roll_num, axis=1)
                p_vel = np.roll(p_vel, roll_num, axis=1)

            p   = np.pi/2
            Vy  = (r_vel * np.sin(theta) * np.sin(p) + (p_vel * np.sin(theta) * np.cos(p) + t_vel * np.cos(theta) * np.sin(p)))
            Vz  = (r_vel * np.cos(p) - p_vel * np.sin(p))
            Vzp = Vy * np.sin(inclin) + Vz * np.cos(inclin)

        
        flux_radius = flux_dens * radius**2 * au_km
        for i in range(len(k_vel)):
            Ij[i] = np.sum(flux_radius * constant2 * np.exp(-(k_vel[i] - Vzp)**2 / b**2))

        max_Ij = np.max(Ij)
        Ij     = Ij / max_Ij


        #----------------------------------------------------------------------------------------------
        #Plot Line Desnsity
        #----------------------------------------------------------------------------------------------

        fig, ax =plt.subplots()
        ax.plot(k_vel,Ij,  label="PEN Inc")
        plt.legend()
        if run_type == "run_stationary":
            plt.title("Stationary Line Density at time " + str(time))
        if run_type == "run_rotate":
            plt.title("Rotational Line Density at time " + str(time))
        plt.xlabel("Velocity [km/s]")
        plt.ylabel("Normalized Flux")
        if dim == '2D':
            plt.savefig("Plots/run_" + str(run_num) + "/" + str(file_name_2D[-1]) + "/" + str(file_name_2D[-1]) + '_' + str(k) + ".png")
        if dim =='3D':
            plt.savefig("Plots/run_" + str(run_num) + "/" + str(file_name_3D[-1]) + "/" + str(file_name_3D[-1]) + '_' + str(k) + ".png")
        plt.close()
        
        #----------------------------------------------------------------------------------------------
        #Polar Plots
        #----------------------------------------------------------------------------------------------
        vel_max = np.max(Vzp)
        vel_min = np.min(Vzp)
        
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(5,5))
        ax.grid(False)
        plt.pcolormesh(data[2],data[1],Vzp)
        ax.set_rmax(0.3)
        ax.set_rticks([0.1, 0.2, 0.3])  # Less radial ticks
        ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
        if run_type == "run_stationary":
            ax.set_title("Projected PEnGUIn Velocities at time " + str(time) + " and Orbit " + str(k))
        if run_type == "run_rotate":
            ax.set_title("Projected PEnGUIn Velocities at time " + str(time) + ",\n Orbit " + str(orbit_num) + ", and Frame " + str(k-frame_start))
        plt.colorbar(label="[km/s]")
        textstr = '\n'.join((
            r'Max Velocity = %.3e'% (vel_max, ),
            r'Min Velocity = %.3e'% (vel_min, )))

        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='gray', alpha=0.3)

        # place a text box in upper left in axes coords
        ax.text(-0.1, -0.09, textstr, transform=ax.transAxes, fontsize=8,
                verticalalignment='bottom', bbox=props)
        plt.savefig("Plots/run_" + str(run_num) + "/" + str(file_name_2D[-3]) + "/polpeng_incvel"  + '_' + str(k) + ".png")
        plt.close()
        
        
        
        dens_max = np.max(density)
        dens_min = np.min(density)

        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.grid(False)
        plt.pcolormesh(data[2],data[1],np.log10(density),vmin=15,vmax=20)
        plt.plot(data[2], radius_1, color='k', linestyle='dotted')
        plt.plot(data[2], radius_2, color='k', linestyle='dotted')
        ax.set_rmax(0.3)
        ax.set_rticks([0.1, 0.2, 0.3])  # Less radial ticks
        ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
        if run_type == "run_stationary":
            ax.set_title("Calculated PEnGUIn Densities at time " + str(time) + " and Orbit " + str(k))
        if run_type == "run_rotate":
            ax.set_title("Calculated PEnGUIn Densities at time " + str(time) + ",\n Orbit " + str(orbit_num) + ", and Frame " + str(k-frame_start))
        plt.colorbar(label="[log(cm-2)]")
        textstr = '\n'.join((
            r'Inner = %.3f' % (inner, ),
            r'Outer = %.3f' % (outer, ),
            r'Max Denisty = %.3e'% (dens_max, ),
            r'Min Denisty = %.3e'% (dens_min, )))

        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='gray', alpha=0.3)

        # place a text box in upper left in axes coords
        ax.text(-0.1, -0.09, textstr, transform=ax.transAxes, fontsize=8,
                verticalalignment='bottom', bbox=props)
        plt.savefig("Plots/run_" + str(run_num) + "/" + str(file_name_2D[-2]) + "/polpeng_dens"  + '_' + str(k) + ".png")
        plt.close()
        
        
        
        #----------------------------------------------------------------------------------------------
        #Multi Plot
        #----------------------------------------------------------------------------------------------
        fig, _axs = plt.subplots(nrows=1, ncols=3, figsize=(18, 5))
        #fig.subplots_adjust(hspace=0.3)
        if run_type == "run_stationary":
            fig.suptitle('Orbit ' +str(k) + ' at Time ' + str(time))
        if run_type == "run_rotate":
            fig.suptitle('Orbit ' +str(orbit_num) + ' at Frame ' + str(k-frame_start) + ' at Time ' + str(time))
        axs = _axs.flatten()

        cset0 = axs[0].plot(k_vel,Ij,  label="PEN Inc")
        cset1 = axs[1].pcolormesh(data[1], data[2], np.log10(np.transpose(density)), vmin=15, vmax=20)
        cset2 = axs[2].pcolormesh(data[1][:-1], data[2][:-1], np.transpose(Vzp), vmin=-150, vmax=150)
        
        if run_type == "run_stationary":
            axs[0].set_title("Stationary Line Density")
        if run_type == "run_rotate":
            axs[0].set_title("Rotational Line Density")
        axs[0].set_xlabel("Velocities [km/s]")
        axs[0].set_ylabel("Normalized Flux")
        axs[1].set_title("Calculated PEnGUIn Densities")
        axs[1].set_xlabel("Radius [au]")
        axs[1].set_ylabel("Theta [rads]")

        axs[2].set_title("PEnGUIn Inclined Velocity")
        axs[2].set_xlabel("Radius [au]")
        axs[2].set_ylabel("Theta [rads]")

        axs[1].plot(radius_1, data[2], color='k', linestyle='dotted')
        axs[1].plot(radius_2, data[2], color='k', linestyle='dotted')

        textstr = '\n'.join((
            r'Inner = %.3f' % (inner, ),
            r'Outer = %.3f' % (outer, ),
            r'Max Denisty = %.3e'% (dens_max, ),
            r'Min Denisty = %.3e'% (dens_min, )))

        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='gray', alpha=0.3)

        # place a text box in upper left in axes coords
        axs[1].text(2.6, -0.14, textstr, transform=ax.transAxes, fontsize=7,
                verticalalignment='bottom', bbox=props)

        textstr = '\n'.join((
            r'Max Velocity = %.3e'% (vel_max, ),
            r'Min Velocity = %.3e'% (vel_min, )))

        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='gray', alpha=0.3)

        # place a text box in upper left in axes coords
        axs[2].text(3.7, -0.12, textstr, transform=ax.transAxes, fontsize=8,
                verticalalignment='bottom', bbox=props)

        fig.colorbar(cset1, ax=axs[1], label="[log(cm-2)]")
        fig.colorbar(cset2, ax=axs[2], label="[km/s]")
        plt.savefig("Plots/run_" + str(run_num) + "/" + str(file_name_2D[-4]) + "/multi_plot"  + '_' + str(k) + ".png")
        plt.close()

        
        #----------------------------------------------------------------------------------------------
        #Save Data
        #----------------------------------------------------------------------------------------------
        save_file_var  = [part_func, flux_dens, opt_depth, Ij, Vzp]
        
        for i in range(len(save_file_name)):
            file_path = "Data/run_" + str(run_num) + "/" + save_file_name[i] + "/data_" + str(k) + "_"
            fname = file_path + save_file_name[i] + ".txt"
            np.savetxt(fname, save_file_var[i])
  
        
    else:
        continue









