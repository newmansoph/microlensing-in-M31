from math import pi, asin
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

###setting up the x-ray binary disc:
theta_E = 1.28*(10**(-5)) #arcsec. for a 0.6 solar mass WD 
d_s = 780 #kpc
d_s_cm = 2.4 * (10**(24)) #cm
R_G = 1474070 #cm
R_G_arcsec = 1.27 * 10**(-13) 
r_E = 6.89*(10**(15)) #cm


###FUNCTIONS:

#if lens inside:
def mag_within_r(r,u0):
    a0 = (4+r*r)**0.5 / r
    a1 = 0
    a2 = -8/(r**3 *(4+r*r)**(3/2))
    a3 = 0
    a4 = -144*(2+2*r*r+ r**4)/(r**5 *(4+r*r)**(7/2))
    mu = a0 + a2*u0*u0/4 + (a4*u0**4)/24
    return mu
    
#if lens outside:
def mag_outside_r(r,u0):
    c0 = (2+u0*u0)/(u0*(4+u0*u0)**0.5)
    c1 = 0
    c2 = (8*(1+u0*u0))/(u0**3 *(4+u0*u0)**(5/2))
    c3 = 0
    c4 = 48*(12+14*u0*u0+6*u0**4+3*u0**6)/(u0**5 *(4+u0**2)**(9/2))
    mu = c0 + (c2*r**2)/4 + (c4*r**4)/24
    return mu 

#for the lens at the source radius:
def mag_at_r(r): #r=xi
    mu = 1/pi * ( 2/r + ((1+r*r)/(r*r)) * ( pi/2 + asin((r*r-1)/(r*r+1)) ))
    return mu

def calc_mags(r_disk,lens_x,lens_y):
    """given the lens trajectory, mags are calculated"""
    
    u0 = []
    for i in np.arange(n):
        u0.append( (lens_x[i]**2 + lens_y[i]**2)**0.5 ) 

    mags = []
    for i in np.arange(n) :

        if lens_x[i] > r_disk :
            mags.append( mag_outside_r(r_disk,u0[i]) )
            #print('outside:',lens_x[i])    
            
        if lens_x[i] < -1*r_disk :
            mags.append( mag_outside_r(r_disk,u0[i]) )
            #print('outside:',lens_x[i]) 
            
        if lens_x[i] < r_disk and lens_x[i] > -1*r_disk :
            mags.append( mag_within_r(r_disk,u0[i]) )
            #rint('inside:',lens_x[i])
           
        if lens_x[i] == r_disk :
            mags.append( mag_at_r(r_disk) )
            #print('at r:',lens_x[i],mags[-1])
            
        if lens_x[i] == -1*r_disk :
            mags.append( mag_at_r(r_disk) )
            #print('at r:',lens_x[i],mags[-1])
    return mags

        
#assuming a WD speed of 200 km/s :
v = 1.38*(10**(-5)) #arcsec/yr

###REGION RADII

#These were obtained with accretion_disk.py: 
r_disk_jwst = 19286*R_G / r_E
r_disk_xmm = 3.64*R_G / r_E 
r_disk_xmm_6Rg = 12.4*R_G / r_E 
r_disk_uvot = 6430*R_G / r_E
r_disk_ztf_g = 8572*R_G / r_E
r_disk_ztf_r = 15001*R_G / r_E
r_disk_ztf_i = 15001*R_G / r_E

##X TRAJECTORIES
x0 = 30000*R_G / r_E
print('x0 = 30,000Rg:',x0, 'Einstein radii')
y0 = 0
r_end = 1.24

###JWST
df1 = pd.read_csv('./data/accretion_outputs/jwst.txt',delimiter=",",usecols=["Total_flux","R_outer"])
flux_left = df1["Total_flux"].tolist() #unmagnified fluxes, one half only
flux_jwst = [value for value in flux_left if value != 0]
del flux_jwst[-1] #removes last value
flux_jwst.extend(flux_left[::-1])  # Extend b with the reversed elements of a
n = len(flux_jwst)
print('n jwst:',n)

lens_x = np.linspace(-x0,x0,n) #units of theta_E
lens_y = []
for i in np.arange(n):
    lens_y.append(y0) 
mags_jwst = calc_mags(r_disk_jwst,lens_x,lens_y)  

jwst_after = []
for i in np.arange(len(mags_jwst)):
    jwst_after.append(mags_jwst[i]*flux_jwst[i])

#jwst_after = interpolate_flux(1,jwst_after)
t_jwst = np.linspace(-x0/v,x0/v,len(jwst_after))
  
    
##SWIFT-UVOT
df2 = pd.read_csv('./data/accretion_outputs/uvot.txt',delimiter=",",usecols=["Total_flux","R_outer"])
flux_left = df2["Total_flux"].tolist() #unmagnified fluxes, one half only
flux_uvot = flux_left.copy()  # Copy the original list
del flux_uvot[-1] #removes last value
flux_uvot.extend(flux_left[::-1])  # Extend b with the reversed elements of a
n = len(flux_uvot)
print('n uvot:',n)

lens_x = np.linspace(-x0,x0,n) #units of theta_E
lens_y = []
for i in np.arange(n):
    lens_y.append(y0) 
mags_uvot = calc_mags(r_disk_uvot,lens_x,lens_y) 

uvot_after = []
for i in np.arange(len(mags_uvot)):
    uvot_after.append(mags_uvot[i]*flux_uvot[i])

t_uvot = np.linspace(-x0/v,x0/v,n)


    
##ztf g-FILTER
df3 = pd.read_csv('./data/accretion_outputs/ztf_g.txt',delimiter=",",usecols=["Total_flux","R_outer"])
flux_left = df3["Total_flux"].tolist() #unmagnified fluxes, one half only
flux_ztf_g = flux_left.copy()  # Copy the original list
del flux_ztf_g[-1] #removes last value
flux_ztf_g.extend(flux_left[::-1])  # Extend b with the reversed elements of a
n = len(flux_ztf_g)
print('n ztf_g:',n)

lens_x = np.linspace(-x0,x0,n) #units of theta_E
lens_y = []
for i in np.arange(n):
    lens_y.append(y0) 
mags_ztf_g = calc_mags(r_disk_ztf_g,lens_x,lens_y) 

ztf_g_after = []
for i in np.arange(len(mags_ztf_g)):
    ztf_g_after.append(mags_ztf_g[i]*flux_ztf_g[i])

t_ztf_g = np.linspace(-x0/v,x0/v,n)    
    

##ztf r-FILTER
df4 = pd.read_csv('./data/accretion_outputs/ztf_r.txt',delimiter=",",usecols=["Total_flux","R_outer"])
flux_left = df4["Total_flux"].tolist() #unmagnified fluxes, one half only
flux_ztf_r = flux_left.copy()  # Copy the original list
del flux_ztf_r[-1] #removes last value
flux_ztf_r.extend(flux_left[::-1])  # Extend b with the reversed elements of a
n = len(flux_ztf_r)
print('n ztf_r:',n)

lens_x = np.linspace(-x0,x0,n) #units of theta_E
lens_y = []
for i in np.arange(n):
    lens_y.append(y0) 
mags_ztf_r = calc_mags(r_disk_ztf_r,lens_x,lens_y)

ztf_r_after = []
for i in np.arange(len(mags_ztf_r)):
    ztf_r_after.append(mags_ztf_r[i]*flux_ztf_r[i])

t_ztf_r = np.linspace(-x0/v,x0/v,n)


##ztf i-FILTER
df5 = pd.read_csv('./data/accretion_outputs/ztf_i.txt',delimiter=",",usecols=["Total_flux","R_outer"])
flux_left = df5["Total_flux"].tolist() #unmagnified fluxes, one half only
flux_ztf_i = flux_left.copy()  # Copy the original list
del flux_ztf_i[-1] #removes last value
flux_ztf_i.extend(flux_left[::-1])  # Extend b with the reversed elements of a
n = len(flux_ztf_i)
print('n ztf_i:',n)

lens_x = np.linspace(-x0,x0,n) #units of theta_E
lens_y = []
for i in np.arange(n):
    lens_y.append(y0) 
mags_ztf_i = calc_mags(r_disk_ztf_i,lens_x,lens_y) 

ztf_i_after = []
for i in np.arange(len(mags_ztf_i)):
    ztf_i_after.append(mags_ztf_i[i]*flux_ztf_i[i])

t_ztf_i = np.linspace(-x0/v,x0/v,n)


###XMM 
df6 = pd.read_csv('./data/accretion_outputs/xmm.txt',delimiter=",",usecols=["Total_flux","R_outer"])

R_outer = df6["R_outer"].tolist()
flux_left = df6["Total_flux"].tolist()
flux_left_new = [value for value in flux_left if value != 0]
    
flux_xmm = flux_left_new.copy()
del flux_xmm[-1] #removes last value
flux_xmm.extend(flux_left_new[::-1])  # Extend b with the reversed elements of a
n = len(flux_xmm)
print('n xmm:',n)

x0 = 30*R_G / r_E
print('x0 = 30Rg:',x0, 'Einstein radii')
y0 = 0
lens_x = np.linspace(-x0,x0,n) #units of theta_E
lens_y = []
for i in np.arange(n):
    lens_y.append(y0) 
mags_xmm = calc_mags(r_disk_xmm,lens_x,lens_y)

xmm_after = []
for i in np.arange(len(mags_xmm)):
    xmm_after.append(mags_xmm[i]*flux_xmm[i])
    
t_xmm = np.linspace(-x0*526000/v,x0*526000/v,len(xmm_after))


###XMM for 6Rg
df7 = pd.read_csv('./data/accretion_outputs/xmm_for_6Rg.txt',delimiter=",",usecols=["Total_flux","R_outer"])

R_outer_6Rg = df7["R_outer"].tolist()
flux_left_6Rg = df7["Total_flux"].tolist()
flux_left_new_6Rg = [value for value in flux_left_6Rg if value != 0]

flux_xmm_6Rg = flux_left_new_6Rg.copy()
del flux_xmm_6Rg[-1] #removes last value
flux_xmm_6Rg.extend(flux_left_new_6Rg[::-1])  # Extend b with the reversed elements of a
n = len(flux_xmm_6Rg)

x0 = 30*R_G / r_E
y0 = 0
lens_x = np.linspace(-x0,x0,n) #units of theta_E
lens_y = []
for i in np.arange(n):
    lens_y.append(y0) 
mags_xmm_6Rg = calc_mags(r_disk_xmm_6Rg,lens_x,lens_y) 

xmm_after_6Rg = []
for i in np.arange(len(mags_xmm_6Rg)):
    xmm_after_6Rg.append(mags_xmm_6Rg[i]*flux_xmm_6Rg[i])
    
t_xmm_6Rg = np.linspace(-x0*526000/v,x0*526000/v,len(xmm_after_6Rg))



###RATIOS

ratio = []
for i in np.arange(len(ztf_r_after)):
    ratio.append(ztf_r_after[i]/jwst_after[i])




###PLOTS 

"""
plt.figure()
plt.plot(t_jwst,ratio,label='ZTF-r/JWST',marker='o')
plt.xlabel('t (years)')
plt.ylabel('Magnified flux')
#plt.yscale('log')
plt.legend()
plt.savefig('./plots/ratio.png',dpi=300,bbox_inches='tight')
"""
    
plt.figure()
plt.plot(t_jwst,jwst_after,label='JWST',marker='o')
plt.plot(t_uvot,uvot_after,label='Swift UVOT',marker='o')
plt.plot(t_ztf_r,ztf_r_after,label='ztf r-filter',marker='o')
plt.xlabel('t (years)')
plt.ylabel('Magnified flux')
plt.yscale('log')
plt.legend()
plt.savefig('./plots/magnified_fluxes.png',dpi=300,bbox_inches='tight')




plt.figure()
plt.plot(t_ztf_g,ztf_g_after,label='g',marker='o')
plt.plot(t_ztf_i,ztf_i_after,label='i',marker='o')
plt.plot(t_ztf_r,ztf_r_after,label='r',marker='o')
plt.xlabel('t (years)')
plt.ylabel('Magnified flux')
plt.yscale('log')
plt.legend()
plt.savefig('./plots/ztf_filters.png',dpi=300,bbox_inches='tight')




plt.figure()
plt.plot(t_xmm,xmm_after,label=r'$r_{isco}=1.24R_g$',marker='o')
plt.plot(t_xmm_6Rg,xmm_after_6Rg,label='$r_{isco}=6R_g$',marker='o')
plt.xlabel('t (minutes)')
plt.ylabel('Magnified flux')
plt.yscale('log')
plt.legend()
plt.savefig('./plots/xmm_peak.png',dpi=300,bbox_inches='tight')




###MAGNIFICATIONS
"""
plt.figure()   
#plt.plot(t_uvot,mags_uvot,label='Swift UVOT',marker='x')
#plt.plot(t_jwst,mags_jwst,label='JWST',marker='x')
plt.plot(t_xmm,mags_xmm,label='XMM',marker='x')
plt.plot(t_xmm_6Rg,mags_xmm_6Rg,label=r'XMM, $6R_g$',marker='x')
plt.xlabel('t (minutes)')
plt.ylabel(r'$\mu$',rotation=0,labelpad=10)
plt.yscale('log')
plt.legend()
plt.savefig('./plots/xmm_mags.png',dpi=300,bbox_inches='tight')

plt.figure()   
plt.plot(t_uvot,mags_uvot,label='Swift UVOT',marker='x')
plt.plot(t_jwst,mags_jwst,label='JWST',marker='x')
plt.xlabel('t (years)')
plt.ylabel(r'$\mu$',rotation=0,labelpad=10)
plt.yscale('log')
plt.legend()
plt.savefig('./plots/witt_mao_mags.png',dpi=300,bbox_inches='tight')
"""





### UNMAGNIFIED FLUXES
plt.figure()
plt.plot(R_outer,flux_left,label=r'$1.24R_g$',marker='x')
plt.plot(R_outer_6Rg,flux_left_6Rg,label=r'$6R_g$',marker='o')
plt.legend()
plt.xlabel('Radius')
plt.ylabel('Unmagnified fluxes (cgs)')
plt.savefig('./plots/unmagnified_fluxes.png',dpi=300,bbox_inches='tight')
