import numpy as np
from scipy.integrate import quad, quadrature, romberg
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["font.size"] = 14

#cgs units 
h = 6.626176*(10**(-27)) #erg second
c = 3*(10**10) #cm/s
k = 1.38*(10**(-16)) #cm2 g s-2 K-1
D = 780 #kpc
D = D * 3.086*(10**21) #cm
G = 6.67*10**(-8) #dyne cm2gm−2.
sigma = 5.6704*10**(-5) # g s-3 K-4
M_Edd = 10**(18) #g/s
M = 10*1.989*(10**33) #g, 10 solar masses

#for a 10 solar mass black hole:
r_E = 6.89*(10**(15)) #cm
R_G = (G*M)/(c*c) #cm
#r_isco = 1.24*R_G #cm, maximal spin for a BH 
r_isco = 6*R_G

f_col = 1.7 

def B(E,T):
    """This is the Planck blackbody spectrum"""
    return (((8*np.pi)/(h*c)**3) * E**3) / (np.exp(E/(k*T))-1)

def integrand(T):
        return (T/T_in)**(-11/3) * B(E,T) / T_in #dT
    

def calc_temp(R):
    return f_col*( (3*G*M_Edd*M)/(8*np.pi*sigma*R**3) * (1 - (r_isco/R)**(1/2) ))**(1/4)

E_start_ev_xmm = 100 #outer radius energy < inner 
E_finish_ev_xmm = 15000 #inner radius energy

E_start_ev_jwst = 0.25
E_finish_ev_jwst = 2.07

E_start_ev_uvot = 1.9
E_finish_ev_uvot = 7.3

E_start_ev_ztf_g = 2.21
E_finish_ev_ztf_g = 3.38

E_start_ev_ztf_r = 1.68
E_finish_ev_ztf_r = 2.26

E_start_ev_ztf_i = 1.78
E_finish_ev_ztf_i = 1.81


#for 1.24Rg I used 25 values
#for 6Rg I set this between 30 and 6 with 16 values
#for JWST I set between 30,000 and 1.24 with 15 values 
#same settings for UVOT and ZTF as JWST 


r_begin = 30  #outer radius to start integrating over
r_end = 6 #inner radius. r_end < r_begin
n_step = 16 #number of radius rings, 300 usually
R_disk = np.linspace(r_begin*R_G,r_end*R_G,n_step) #cm

E_start_ev = E_start_ev_xmm
E_finish_ev = E_finish_ev_xmm

ev_to_erg = 1.60218*10**(-12)
E_start_erg = E_start_ev *ev_to_erg #converting eV to ergs
E_finish_erg = E_finish_ev * ev_to_erg
n_int = 16 #number of energy steps for the integration
E_values = np.linspace(E_start_erg,E_finish_erg,n_int) #for optical band 
E_values_ev = np.linspace(E_start_ev,E_finish_ev,n_int)


total_flux = []
integrand_to_plot = []

f = open('./data/accretion_outputs/xmm_for_6Rg.txt','w')  
f.close()
file = open('./data/accretion_outputs/xmm_for_6Rg.txt','a+') #opens file in append and read mode
file.write('R_outer,R_inner,Total_flux')

"""for each small radii range in the energy band"""
for i in np.arange(len(R_disk)-1):

    R_outer = R_disk[i]
    R_inner = R_disk[i+1]

    #temps that we are integrating over:
    T_values = []
    R_values = np.linspace(R_outer,R_inner,n_int)

    for w in R_values: 
        T_values.append( calc_temp(w) )

    print(i)

    flux = []
    
    E_to_integrate = []
    #performing the integral:
    for p in np.arange(len(T_values)-1):
        T_in = T_values[p+1]
        T_out = T_values[p]
        T_to_integrate = np.arange(T_out,T_in,n_int)
        
        E = (E_values[p]+E_values[p+1])/2   
        E_to_integrate.append(E)
        
        r_in = R_values[p+1]
        before_integrand = (8*np.pi*(r_in)**2)/(3*D*D)
        #integral = quad(integrand,T_out,T_in)[0]
        
        integrand = []
        for point in T_to_integrate:
            integrand.append( (point/T_in)**(-11/3) * (1/T_in) * B(E,point) )
        
        integral = np.trapz(integrand, T_to_integrate)
        
        ##if T_in < T_out :
          ##  integral = np.trapz(integrand,T_in,T_out)[0]
     
        flux.append(before_integrand*integral) #in each smaller slice

    total_flux.append( np.trapz(x=E_values_ev[0:-1],y=flux) ) #in each big slice 
    file.write(f'\n{R_outer/R_G},{R_inner/R_G},{total_flux[-1]}')


"""
I recommend plotting this when running the code to check everything looks right with the energy/radii steps you have chosen
plt.figure()
plt.plot(R_disk[0:n_step-1]/R_G, total_flux,ls='None',marker='o',markersize=2) #for 100 R values there is 99 flux values
plt.xlabel('r')
plt.ylabel('Flux')
#plt.yscale('log')
plt.savefig('./plots/accretion_disk/flux_vs_r.png',dpi=300, bbox_inches="tight")   
"""



###making the heatmap:



R_disk = np.linspace(r_begin,r_end,n_step)

num_extra_points = 2

# Interpolate values between original data points
interpolated_flux = []
interpolated_radii = []
for i in range(len(total_flux) - 1):
    flux_start = total_flux[i]
    flux_end = total_flux[i + 1]
    radii_start = R_disk[i]
    radii_end = R_disk[i + 1]
    
    # Linear interpolation
    flux_values = np.linspace(flux_start, flux_end, num_extra_points + 2)
    radii_values = np.linspace(radii_start, radii_end, num_extra_points + 2)
    
    # Add interpolated values to the lists
    interpolated_flux.extend(flux_values)
    interpolated_radii.extend(radii_values)

fig = plt.figure()
ax = plt.subplot(projection="polar")

azm = np.linspace(0, 2 * np.pi, len(interpolated_radii))
r, th = np.meshgrid( np.resize(interpolated_radii,len(azm)) , azm)

###where does max flux occur?
index_max = interpolated_flux.index(max(interpolated_flux))
r_max = interpolated_radii[index_max]
print(r_max)

#total_flux = np.array([flux for flux in total_flux] + [total_flux[-1]]) #repeating the last value to allow plotting

z= np.resize(interpolated_flux,[len(r),len(th)])

r_labels=[r_end,r_max,r_begin]
rticks = ax.set_rticks(r_labels)

#to deal with XMM label not showing:
if E_start_ev == E_start_ev_xmm :
    rlabels = ax.get_ymajorticklabels()
    for label in rlabels:
        label.set_color('white')

ax.set_xticks([]) #removes theta labels




c = plt.pcolormesh(th,r,z, shading='nearest',cmap="gist_heat")
clb = plt.colorbar(c,pad=0.2)
clb.set_label("Flux (erg/cm$^2$/s)",labelpad=10) #rotation=270

plt.plot(azm, r, color='k', ls='none') 
plt.suptitle(f'XMM')


plt.savefig('./plots/accretion_disk/xmm_6rg.png',dpi=300,bbox_inches='tight')

print('max flux:',interpolated_flux[index_max])
