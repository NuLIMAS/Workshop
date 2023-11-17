import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import splev, splrep

# Function to calculate period-averaged pressure
def calculate_period_averaged_pressure(T, t_instantaneous, p_instantaneous):
    num_cycles = int(np.floor(t_instantaneous[-1] / T))
    mean_pressures = []
    time_instances = []
    
    for i in range(num_cycles):
        cycle_start_index = np.where(t_instantaneous >= i * T *0.95)[0][0]
        cycle_end_index = np.where(t_instantaneous <= (i + 1) * T*1.05)[0][-1]
        model = splrep(t_instantaneous[cycle_start_index:cycle_end_index+1], p_instantaneous[cycle_start_index:cycle_end_index+1])
        presure_avg = np.mean(splev(np.linspace(i * T,  (i+1) * T,25), model ))
        mean_pressures.append(presure_avg)
        time_instances.append((i+1)*T)

    
    
    return mean_pressures, time_instances

# Constants
T = 1.6
H = 0.15
depths = [0.05, 0.1, 0.15, 0.18]
alphainit = 0.5 # Set initial porosity
rhos = 2656
rhof = 1000
g = 9.81

#Size settings
fontsizelegend =20
fontsizelabel =20
ticksize =20
titlesize =16

# File paths
file_path_num = f'postProcessing/Probes/0/p'


# Load data
data_num = np.loadtxt(file_path_num)


t_instantaneous_num = data_num[:, 0]
p_instantaneous_num = [data_num[:, i] / 1000 for i in [1, 2, 3, 4]]

# Create a plot for pore pressure
plt.figure(figsize=(10, 6))
colors = ['blue', 'green', 'red', 'purple']
max_mean_pressures_num = [0]



for i, (p_num) in enumerate( p_instantaneous_num):
    mean_pressures_num, time_instance = calculate_period_averaged_pressure(T, t_instantaneous_num, p_num)
    #print(time_instance )
    #plt.plot(np.linspace(0, t_instantaneous_num[-1], len(mean_pressures_num)), mean_pressures_num, linestyle='-', color=colors[i], label=f'Z = {depths[i]} m')
    plt.plot(time_instance, mean_pressures_num, linestyle='-', color=colors[i], label=f'Z = {depths[i]} m')
    rho = alphainit * rhos + (1 - alphainit) * rhof
    

plt.xlabel('Time (s)',fontsize=fontsizelabel)
plt.ylabel('Pressure (kPa)',fontsize=fontsizelabel)
#plt.title(f'Pore Pressure for T = {T} s and H = {H:.2f} m',fontweight='bold',fontsize=titlesize)
legend=plt.legend()
for text in legend.get_texts():
    #text.set_fontweight('bold')
    text.set_fontsize(fontsizelegend)
#plt.ylim(0, 7)
#plt.xlim(0, 4000)
plt.xticks(fontsize=ticksize,)
plt.yticks(fontsize=ticksize)
plt.savefig(f'pore_pressure_T_{T}_H_{H:.2f}.png')
plt.show()



