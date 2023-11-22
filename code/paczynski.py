"""this code plots the paczynski curves using the results of sim.py"""

import pandas
from matplotlib import pyplot as plt
import numpy as np
import random
from scipy.optimize import curve_fit

import style
plt.style.use(style.custom_style)

d_m31 = 780 #kpc
d_m31_cm = d_m31 * 3.086e+21 #cm
    
def min_flux(exp_time):
    """The detection sensitivity of Swift is 2 × 10−14 erg/cm2/s in 10^4"""
    return( (10**4)*2*(10**(-14)) / exp_time )

#average swift exposure time 
xamin = pandas.read_csv('/home/sophienewman/microlensing/my_paper/data/xamin.csv',usecols=["xrt_exposure","start_time"])

swift_exp = xamin["xrt_exposure"].tolist()
swift_exp_av = 2270
swift_times = xamin["start_time"] #in UTC format

#flux_errors = pandas.read_csv('/home/sophienewman/microlensing/my_paper/data/flux_errors.csv',usecols=["Error on luminosity"])
swift_exp = xamin["xrt_exposure"].tolist()

sw_flux_min = min_flux(swift_exp_av) 
print('min sw flux for 2270s:',sw_flux_min)
sw_lum_min = sw_flux_min*4*np.pi*d_m31_cm*d_m31_cm
sw_lum_min = 5.4e+36

obs_seps_mean = 1 #days
obs_seps_mean = obs_seps_mean / 365 #years 
print(obs_seps_mean,'years')


def indices(lst, item):
    return [i for i, x in enumerate(lst) if x == item]
    
def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs
    
def closest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx




def paczynski_curve(t, u0, t_cross):
    u_t = np.sqrt(u0 * u0 + (t / t_cross) ** 2)
    return (u_t * u_t + 2) / (u_t * np.sqrt(u_t * u_t + 4))



### SRC LUMINOSITIES 
""" ###this uses the 164 flux values I got from my merged image
chandra_data = pandas.read_csv('/home/sln/mphys_lensing/filtered.csv',usecols=["RA","DEC","COUNTS ","COUNT RATE","FLUX (ERG/CM/CM/S)"])
ch_flux = chandra_data["FLUX (ERG/CM/CM/S)"].tolist()
ch_rate = chandra_data["COUNT RATE"]
ch_counts = chandra_data["COUNTS "]
ch_lum = []
for f in ch_flux:   
    ch_lum.append( 4*np.pi*d_m31_cm*d_m31_cm*f ) #erg/s
"""  

chandra_data = pandas.read_csv('/home/sophienewman/microlensing/my_paper/data/m31_catalogue_within_fov.txt',usecols=["Summed flux (0.5-8.0 keV)","net cts","Tot Exp (ks)"])
ch_flux = chandra_data["Summed flux (0.5-8.0 keV)"].tolist()
ch_counts = chandra_data["net cts"].tolist()
ch_exp = chandra_data["Tot Exp (ks)"].tolist() #units of ks

# Remove the dodgy values from the lists
indices_to_remove = np.where(ch_flux == -9.99E+00)[0]
for i in indices_to_remove:
    ch_flux.pop(i)
    ch_counts.pop(i)
    
    
    
#- convert Chandra count rate to Swift count rate assuming flat spectrum

r2_ch = ch_counts[0] / (ch_exp[0]*1000) #0.0011357142857142857
r2_sw = 5.319E-04  #found using WebPimms and r2_ch

#- linearly scale count rate with flux 
f2_ch = ch_flux[0] #from Chandra data, 2.24e-14
f1_sw = 4.783E-14 #from WebPimms using f2_ch

r1_sw = (f1_sw/f2_ch) * r2_sw #gives 0.001 ish

    
ch_lum = []
for f in ch_flux:   
    ch_lum.append( 4*np.pi*d_m31_cm*d_m31_cm*f ) #erg/s
    
    

n_sims = 94
outputlist = np.arange(1,n_sims-1)
outputlist = outputlist.tolist()

    

###READING OUTPUT FILES 

for output_file_number in outputlist :

    data1 = pandas.read_csv(f'/home/sophienewman/microlensing/my_paper/sim_results/outputs{output_file_number}.txt',sep=',',usecols=["lens","u0","A","delta_t"])
    
    n_lenses = data1["lens"]
    n_events = list(dict.fromkeys(n_lenses)) #removing duplicates 
    
    lenses_ran = 3000000
    
    ### PROBABILITY
    
    prob = len(n_events) / lenses_ran
    
    n_above = 0
    
    for event in np.arange(len(n_events)): #for each event 
    
        n = n_events[event] #one of my lensing events 
    
        evt = data1.loc[data1["lens"]==n]
        u0_list = evt["u0"].tolist()
        if u0_list[0] != 'u0': 
            u0_values = [n for n in u0_list if n < 1] #list of u0 only less than 1
            
            if len(u0_values) > 1:
                
                #for the data points:
                u0 = min(u0_values)
                u0_idx = u0_list.index(u0)
            
                t_cross = evt["delta_t"].tolist()
                t_cross = abs(t_cross[u0_idx])
                t_cross = np.sqrt(780/10) * t_cross #for M31 lenses
                t_cross_rounded = round(t_cross,2)
                
                #Paczynski curve using u0 and t_crossing:
                t = np.linspace(-0.5*t_cross,0.5*t_cross,len(u0_values))
                A = []
                for i in np.arange(len(t)):
                    u_t = np.sqrt(u0_values[i]*u0_values[i] + (t[i]/t_cross)**2)
                    A.append( (u_t*u_t + 2) / (u_t*np.sqrt(u_t*u_t + 4)) )
                             
                # Use curve_fit to fit the data
                popt, pcov = curve_fit(paczynski_curve, t, A,p0=(1.0, 1.0))

                # popt contains the optimized values for the parameters (u0, t_cross)
                u0_fit, t_cross_fit = popt
                t_fit = np.linspace(-0.5*t_cross,0.5*t_cross,1000)
                A_fit = paczynski_curve(t_fit, u0_fit, t_cross_fit)       
            
                #picking a random source for this event to correspond to:
                index = random.choice(np.arange(len(ch_flux)))
                
                luminosity = []
                flux = []
                counts = []
                rate = []
                for mag in A_fit:
                    luminosity.append( ch_lum[index] * (mag) )
                    flux.append( ch_flux[index] * (mag) )
                    r2 = r1_sw * ch_flux[index]/f1_sw * (mag)  #count rate
                    rate.append(r2)
                    counts.append(r2 * swift_exp_av)
                    
            
                num_swift_line = 100
                swift_line_x = np.linspace(-0.5*t_cross,0.5*t_cross,num=num_swift_line)
                swift_line_y = []
                for n in np.arange(num_swift_line):
                    swift_line_y.append( sw_flux_min ) #plotting Swift's sensitivity
                
                
                #FINDING CROSSING TIMES 
                
                idx1 = 0
                idx2 = len(flux) - 1
                
                # Find the first index above sw_flux_min and last index before it falls below
                for i in range(len(flux)):
                    if flux[i] > sw_flux_min:
                        if idx1 == 0:
                            idx1 = i
                            idx2 = len(flux) - idx1 - 1 #because the Paczynski curve is symmetric:
                            break  # Exit the loop once idx1 and idx2 are found
                           
              
                plt.figure()
                
                
                #if max(flux) > sw_flux_min:
                
                #sw_cross = abs(t_fit[idx2] - t_fit[idx1])
                
                #print('sw_cross:',sw_cross,'swift sep:',obs_seps_mean)
                #num_sampled = int(sw_cross/obs_seps_mean)
                #print('num_sampled:',num_sampled)
                
                # Calculate the sampling rate if num_sampled is non-zero
                #indices = np.linspace(idx1, idx2, num_sampled, dtype=int)
            
                errors = []
                for sample in np.arange(len(flux)) :
                    errors.append( flux[sample] * counts[sample]**(-1/2)  ) #/ (swift_exp_av)
                
        
                
                plt.scatter(t_fit,flux,s=5,zorder=4,alpha=0.6) 
                plt.errorbar(t_fit,flux,yerr=errors,ls='None',zorder=2,alpha=0.3)
                    
                #plt.plot(t_fit,flux,ls='--',zorder=3) #Paczynski curve, color='dodgerblue'
            
                plt.plot(swift_line_x,swift_line_y,ls='--',zorder=1, color='black') 
                
                
                plt.xlabel('t (years)')
                plt.ylabel(f'Flux (erg/s/cm$^2$)')
                
                # Find the mean flux value to center the y-axis
                mean_flux = np.mean(flux)
                y_axis_range = max(flux) - min(flux)
                half_range = y_axis_range / 2
                plt.ylim(bottom= min(flux) - 0.5*half_range, top= max(flux) + 0.5*half_range)

                plt.suptitle(f'A = {round(max(A),2)}, $u_0$ = {round(u0,2)}, $t_E$ = {t_cross_rounded} yrs', y=0.95)
                plt.savefig(f'/home/sophienewman/microlensing/my_paper/sim_analysis/Swift/mags{output_file_number}_{event+1}.png',dpi=300,bbox_inches='tight')

                if max(flux) < swift_line_y[0]:
                    print(output_file_number,event+1)
