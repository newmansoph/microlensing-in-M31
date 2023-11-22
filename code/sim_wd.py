#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
from matplotlib import cm
import random
import math
from statistics import variance
from astropy import units as u
from astropy.coordinates import SkyCoord


plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["font.size"] = 16


# In[3]:

n_lenses = 3000000
n_d = 1000

#distances:

m31 = 780 #kpc
d0 = 10 #kpc
d_halo_mw = 100 #kpc
n_step = (d_halo_mw-d0) / n_d
d_mw = np.arange(d0,d_halo_mw,n_step)

d0_m31 = m31 - 30
d1_m31 = m31 - 5
n_step = (d1_m31 - d0_m31) / n_d
d_m31 = np.arange(d0_m31,d1_m31,n_step)



#for 2D Gaussian (values obtained in wd_dist code):

v_mu = 194.65 # km/s
M_mu = 0.58 #Msolar
mu = np.array([M_mu,v_mu]) # 2 dimensional mean vector

var_M = 0.00271
var_v = 3820
cov_vM = -0.401

cov = np.array([[var_M,cov_vM],[cov_vM,var_v]])  # 2x2 covariance matrix 


### sampling 100 (n_lenses) velocities and masses from this Gaussian to get a distance probability density with 100 distances:

sample = np.random.multivariate_normal(mu,cov,n_lenses)     
mass = sample[:,0]
v = sample[:,1]


v = v.tolist()

#plt.figure()
#plt.scatter(mass,v,alpha=0.1)
#plt.xlabel('Mass/$M_{\odot}$')
#plt.ylabel('Velocity (km/s)')
#plt.savefig('/home/sln/summer_2023/plots/wd_dist_check.png',dpi=300,bbox_inches='tight')



###using density profile of WDs as seen in 
###Ruiter 2007, Zinn 1985, Morrison 1996, Morrison & Sarajedini  1996, Siegel 2002:

dNdt = 1.03*(10**7) #/Myr
dt =  11000 #Myrs, age of halo 
n_mw = [] #black holes/area 
n_m31 = []
a0 = 3.5 #kpc 
count = 0
for i in d_mw:
    n_mw.append( (dNdt*dt) / (1+i/a0)**3.5 )
#for i in d_m31 :
#    n_m31.append( (dNdt*dt) * (1+i/a0)**3.5 )
    

#plt.figure()
#plt.plot(d_mw,n_mw,color='lightseagreen',marker='o',ls='None')
#plt.ylabel('Number of lenses per area')
#plt.xlabel('Distance (kpc)')
#plt.savefig('/home/sln/summer_2023/plots/numberoflenses.png',dpi=300,bbox_inches='tight')
#plt.close()

   
    
#probability densities to sample d from that uniform distance distribution:
prob_mw = []
prob_m31 = []
for i in n_mw:
    prob_mw.append( i/np.sum(n_mw) )
#for i in n_m31:
#    prob_m31.append( i/np.sum(n_m31) )


#my source regions:
ellipses = [ "ellipse(0:41:48.4196,+41:01:18.156,0.49951,0.44312,152.61796)",
"ellipse(0:42:32.9631,+41:03:29.655,0.63551,0.45818,18.021382)",
"ellipse(0:42:29.3457,+41:04:35.862,0.33641,0.25148,40.812419)",
"ellipse(0:42:07.5183,+41:04:39.479,0.45254,0.41987,173.66944)",
"ellipse(0:41:43.1545,+41:05:07.549,0.28112,0.22828,19.651605)",
"ellipse(0:42:38.0202,+41:05:27.030,0.49749,0.40932,60.957579)",
"ellipse(0:42:10.9354,+41:06:47.341,0.97321,0.44814,19.694724)",
"ellipse(0:42:33.4379,+41:06:50.741,0.25664,0.24596,122.00181)",
"ellipse(0:42:41.0191,+41:07:01.470,0.32358,0.21118,36.852328)",
"ellipse(0:42:23.3189,+41:07:37.502,0.23761,0.20212,172.08543)",
"ellipse(0:42:04.2494,+41:09:31.063,0.26016,0.23149,126.69863)",
"ellipse(0:42:28.2070,+41:10:00.329,0.20152,0.18064,15.26491) ",
"ellipse(0:43:03.1259,+41:10:16.114,0.32923,0.30481,4.280112) ",
"ellipse(0:42:07.3854,+41:10:28.524,0.40728,0.34057,144.39024)",
"ellipse(0:42:51.1594,+41:10:29.686,0.31387,0.19954,36.557757)",
"ellipse(0:42:40.7377,+41:10:34.351,0.15351,0.11248,123.52772)",
"ellipse(0:43:32.3378,+41:10:39.948,0.25611,0.21446,164.96009)",
"ellipse(0:42:11.7767,+41:10:48.551,0.19972,0.17505,5.824894) ",
"ellipse(0:42:57.9156,+41:11:04.755,0.12830,0.10055,24.858945)",
"ellipse(0:41:57.2415,+41:11:09.443,0.30464,0.22535,28.013078)",
"ellipse(0:41:44.7445,+41:11:11.335,0.39236,0.34193,24.527051)",
"ellipse(0:42:47.7816,+41:11:13.813,0.12563,0.09602,135.15663)",
"ellipse(0:42:05.4690,+41:11:37.120,0.28050,0.20717,57.849593)",
"ellipse(0:42:44.8182,+41:11:38.344,0.10324,0.07505,6.749191) ",
"ellipse(0:42:47.2229,+41:11:58.108,0.09736,0.06779,10.60423) ",
"ellipse(0:43:01.9053,+41:12:01.721,0.14037,0.11994,21.890562)",
"ellipse(0:42:44.3855,+41:11:58.316,0.06788,0.05209,9.508051) ",
"ellipse(0:41:50.3593,+41:12:12.811,0.29758,0.23887,31.35619) ",
"ellipse(0:42:28.2858,+41:12:22.785,0.09908,0.09769,9.627881) ",
"ellipse(0:42:18.3185,+41:12:23.428,0.15129,0.12704,176.73227)",
"ellipse(0:43:43.7537,+41:12:34.038,0.51657,0.42703,73.529112)",
"ellipse(0:42:15.1424,+41:12:34.043,0.16789,0.13241,176.69768)",
"ellipse(0:42:59.4390,+41:12:42.043,0.18868,0.13187,11.641376)",
"ellipse(0:43:08.5591,+41:12:49.028,0.17155,0.10099,178.97792)",
"ellipse(0:42:10.9374,+41:12:48.838,0.20170,0.18013,174.47661)",
"ellipse(0:42:52.4840,+41:13:01.696,0.18376,0.13001,126.72455)",
"ellipse(0:42:32.1693,+41:13:13.657,0.14520,0.08760,30.12663) ",
"ellipse(0:43:34.2628,+41:13:22.690,0.23254,0.19655,161.25312)",
"ellipse(0:42:52.6888,+41:13:28.727,0.12093,0.10574,56.646144)",
"ellipse(0:42:40.6037,+41:13:28.151,0.10809,0.09911,28.74817) ",
"ellipse(0:42:22.3940,+41:13:33.373,0.11466,0.09301,177.76936)",
"ellipse(0:42:25.1415,+41:13:40.636,0.10016,0.08601,8.067491) ",
"ellipse(0:42:36.6730,+41:13:49.285,0.10766,0.07355,91.433772)",
"ellipse(0:42:18.6189,+41:14:01.306,0.12903,0.10562,11.538209)",
"ellipse(0:42:47.1495,+41:14:07.882,0.08591,0.06426,90.981898)",
"ellipse(0:42:23.1246,+41:14:07.767,0.08576,0.07917,72.756794)",
"ellipse(0:42:21.5606,+41:14:19.605,0.10233,0.07615,1.737796) ",
"ellipse(0:42:39.4658,+41:14:27.817,0.07933,0.04248,141.0984) ",
"ellipse(0:43:37.2636,+41:14:42.568,0.25613,0.21200,161.76303)",
"ellipse(0:42:42.2682,+41:14:44.759,0.09547,0.04979,90.071057)",
"ellipse(0:43:10.5706,+41:14:51.217,0.10769,0.07503,102.39064)",
"ellipse(0:42:10.2734,+41:15:10.075,0.16894,0.15745,138.75802)",
"ellipse(0:42:43.5946,+41:15:14.056,0.16486,0.15914,30.900247)",
"ellipse(0:43:02.9959,+41:15:24.384,0.15416,0.07180,128.55959)",
"ellipse(0:42:48.5437,+41:15:21.197,0.10250,0.03384,89.66774) ",
"ellipse(0:42:45.0470,+41:15:23.263,0.11381,0.10701,142.96682)",
"ellipse(0:42:41.8093,+41:15:27.272,0.22062,0.08316,123.61926)",
"ellipse(0:42:58.3178,+41:15:29.635,0.10234,0.05385,86.696731)",
"ellipse(0:42:22.9665,+41:15:35.051,0.09854,0.07241,5.168027) ",
"ellipse(0:42:52.4095,+41:15:38.801,0.07479,0.05546,174.2916) ",
"ellipse(0:42:40.0030,+41:15:46.737,0.10979,0.05776,179.16248)",
"ellipse(0:42:42.5655,+41:15:54.509,0.08920,0.06428,177.67956)",
"ellipse(0:43:04.2541,+41:16:01.955,0.07715,0.06305,89.620347)",
"ellipse(0:42:54.8445,+41:16:02.438,0.07916,0.06890,13.831344)",
"ellipse(0:42:49.2377,+41:16:02.237,0.07599,0.06250,93.404427)",
"ellipse(0:42:38.5309,+41:16:02.713,0.10979,0.08029,2.238015) ",
"ellipse(0:42:21.4692,+41:16:01.436,0.11747,0.09987,10.86587) ",
"ellipse(0:42:59.7887,+41:16:06.275,0.11633,0.07112,12.148885)",
"ellipse(0:42:44.3219,+41:16:08.806,0.24308,0.14315,118.09546)",
"ellipse(0:42:33.8806,+41:16:18.901,0.09098,0.04642,91.56746) ",
"ellipse(0:42:31.0940,+41:16:22.006,0.05735,0.05642,55.9587)  ",
"ellipse(0:42:47.1494,+41:16:29.379,0.07698,0.03507,89.438665)",
"ellipse(0:42:43.7341,+41:16:32.469,0.12227,0.06677,51.021969)",
"ellipse(0:43:14.2534,+41:16:50.423,0.39052,0.26120,152.17237)",
"ellipse(0:42:52.3845,+41:16:49.539,0.06926,0.06550,56.449559)",
"ellipse(0:43:53.5518,+41:16:54.971,0.40151,0.37057,156.24993)",
"ellipse(0:42:30.3561,+41:16:53.515,0.11132,0.08724,3.718902) ",
"ellipse(0:42:42.5995,+41:16:58.276,0.12083,0.09189,93.551849)",
"ellipse(0:43:24.8283,+41:17:27.220,0.20221,0.16464,84.70204) ",
"ellipse(0:42:52.3477,+41:17:39.824,0.28169,0.24545,84.603588)",
"ellipse(0:42:44.9338,+41:17:40.634,0.11633,0.06731,1.586965) ",
"ellipse(0:42:33.1910,+41:17:42.288,0.09915,0.04446,88.868331)",
"ellipse(0:43:21.0792,+41:17:49.892,0.16841,0.14328,96.69171) ",
"ellipse(0:43:03.8915,+41:18:04.581,0.11485,0.10213,94.376383)",
"ellipse(0:43:11.1909,+41:18:08.973,0.32242,0.21684,3.431009) ",
"ellipse(0:42:27.8816,+41:18:20.636,0.27962,0.22032,53.481747)",
"ellipse(0:43:13.2507,+41:18:13.093,0.12138,0.08962,121.17824)",
"ellipse(0:42:49.2421,+41:18:16.260,0.04493,0.04451,178.3898) ",
"ellipse(0:43:27.8731,+41:18:28.633,0.19939,0.17379,152.16716)",
"ellipse(0:42:55.3360,+41:18:35.804,0.11032,0.05865,11.667746)",
"ellipse(0:43:16.1568,+41:18:40.838,0.57445,0.23258,98.957916)",
"ellipse(0:42:56.9430,+41:18:43.835,0.15595,0.12924,74.438818)",
"ellipse(0:42:40.1703,+41:18:44.492,0.10734,0.08157,91.302414)",
"ellipse(0:42:52.4496,+41:18:54.897,0.09004,0.08309,125.95215)",
"ellipse(0:43:09.8980,+41:19:01.779,0.13935,0.10992,96.341935)",
"ellipse(0:43:26.2771,+41:19:08.545,0.21558,0.19948,48.312375)",
"ellipse(0:42:46.9197,+41:19:11.212,0.15553,0.13350,4.319247) ",
"ellipse(0:42:31.3237,+41:19:10.338,0.28803,0.18471,168.27655)",
"ellipse(0:43:06.7026,+41:19:12.895,0.23144,0.17807,79.986661)",
"ellipse(0:42:26.0921,+41:19:15.170,0.11407,0.09893,123.78028)",
"ellipse(0:42:59.6922,+41:19:19.293,0.08184,0.08064,96.81084) ",
"ellipse(0:41:41.5370,+41:19:18.059,0.51477,0.41516,114.15399)",
"ellipse(0:43:14.6441,+41:19:29.685,0.27563,0.23212,71.319488)",
"ellipse(0:42:31.2741,+41:19:38.896,0.12370,0.09950,156.85532)",
"ellipse(0:41:44.4073,+41:19:50.497,0.27552,0.20736,49.466029)",
"ellipse(0:42:35.2503,+41:20:05.786,0.10925,0.10186,82.716888)",
"ellipse(0:43:07.4495,+41:20:19.086,0.18602,0.15824,96.40814) ",
"ellipse(0:42:15.6158,+41:20:32.952,0.16279,0.15516,121.86865)",
"ellipse(0:42:54.1103,+41:20:35.353,0.17015,0.11829,100.27228)",
"ellipse(0:43:02.9603,+41:20:42.148,0.17401,0.15579,115.85928)",
"ellipse(0:42:27.6484,+41:20:47.668,0.23715,0.19702,86.023917)",
"ellipse(0:42:08.9845,+41:20:48.849,0.37979,0.29804,12.711434)",
"ellipse(0:42:41.6307,+41:21:05.533,0.15285,0.11924,120.06856)",
"ellipse(0:43:50.8476,+41:21:16.627,0.62593,0.43607,110.15861)",
"ellipse(0:43:39.3280,+41:21:15.188,0.56229,0.39134,163.70597)",
"ellipse(0:42:46.8632,+41:21:18.793,0.11174,0.09897,116.18792)",
"ellipse(0:41:50.8275,+41:21:12.867,0.36531,0.30058,113.66668)",
"ellipse(0:43:03.2647,+41:21:21.412,0.17937,0.15395,74.15828) ",
"ellipse(0:42:34.1903,+41:21:49.166,0.19543,0.13228,120.6845) ",
"ellipse(0:42:19.7464,+41:21:56.278,0.37817,0.24878,63.645289)",
"ellipse(0:43:56.3297,+41:22:02.366,0.45863,0.42071,147.38215)",
"ellipse(0:42:52.3927,+41:22:04.741,0.34354,0.26410,161.26535)",
"ellipse(0:42:49.0175,+41:24:06.933,0.26085,0.21936,133.84646)",
"ellipse(0:43:44.5056,+41:24:08.972,0.59644,0.47641,92.852801)",
"ellipse(0:42:45.9611,+41:24:32.395,0.17817,0.16988,68.931203)",
"ellipse(0:42:43.6833,+41:25:18.609,0.26258,0.21008,135.73625)",
"ellipse(0:42:48.4679,+41:25:22.078,0.27289,0.22380,143.95539)",
"ellipse(0:42:26.1751,+41:25:51.669,0.29086,0.26268,133.90427)",
"ellipse(0:42:55.3258,+41:25:56.841,0.30052,0.24900,150.0928) ",
"ellipse(0:41:45.6623,+41:26:18.694,0.59732,0.55344,63.498659)",
"ellipse(0:42:51.1471,+41:26:38.466,0.29757,0.22689,34.968573)",
"ellipse(0:42:20.4809,+41:26:39.959,0.54054,0.41916,175.04693)",
"ellipse(0:42:41.9717,+41:28:31.262,0.27248,0.25956,137.30395)",
"ellipse(0:42:29.4325,+41:29:02.649,0.43407,0.40823,159.76172)",
"ellipse(0:42:18.7973,+41:29:25.364,0.76422,0.51422,167.83686)",
"ellipse(0:42:52.0244,+41:31:08.312,0.52275,0.45984,144.50144)",
"ellipse(0:42:15.6351,+41:31:08.365,0.38996,0.35209,75.548109)",
"ellipse(0:42:34.5471,+41:32:50.265,0.65536,0.57562,138.71233)",
"ellipse(0:43:12.5936,+41:02:07.797,0.17663,0.16069,180)      ",
"ellipse(0:41:34.7242,+41:03:05.715,0.30808,0.23702,134.55241)",
"ellipse(0:41:58.4062,+41:03:29.880,0.16548,0.13747,148.57484)",
"ellipse(0:41:41.3259,+41:03:30.337,0.50821,0.44255,34.398225)",
"ellipse(0:42:15.9250,+41:06:03.936,0.45704,0.37109,17.527295)",
"ellipse(0:41:17.2776,+41:06:34.346,0.27691,0.24635,71.899524)",
"ellipse(0:42:19.0514,+41:08:04.924,0.25029,0.21341,20.223116)",
"ellipse(0:41:55.9091,+41:10:15.353,0.27597,0.25494,76.763718)",
"ellipse(0:41:40.6516,+41:10:36.701,0.16050,0.12597,128.7356) ",
"ellipse(0:43:07.3833,+41:10:58.282,0.26972,0.24600,155.52845)",
"ellipse(0:43:57.8619,+41:13:24.921,0.28699,0.26650,156.84605)",
"ellipse(0:42:01.2211,+41:14:01.307,0.58204,0.31133,36.668594)",
"ellipse(0:44:04.3265,+41:17:18.461,0.22096,0.12859,45.716048)",
"ellipse(0:41:36.2051,+41:17:49.178,0.64575,0.33666,63.538203)",
"ellipse(0:44:02.8135,+41:19:52.459,0.59175,0.43594,80.910328)",
"ellipse(0:41:21.9710,+41:21:40.082,0.56588,0.27398,160.03496)",
"ellipse(0:43:00.4130,+41:22:27.907,0.29751,0.22654,44.173125)",
"ellipse(0:44:00.3357,+41:22:56.852,0.59651,0.42406,28.76542) ",
"ellipse(0:42:59.2284,+41:27:11.363,0.42921,0.28691,2.920601) ",
"ellipse(0:42:46.3550,+41:27:26.850,0.26974,0.22875,76.489814)",
"ellipse(0:42:58.3767,+41:28:37.949,0.29321,0.24288,156.36545)",
"ellipse(0:42:16.0320,+41:28:45.630,0.28467,0.22568,62.968404)",
"ellipse(0:42:53.7856,+41:29:45.630,0.41310,0.36622,117.76758)",
"ellipse(0:42:54.8827,+41:32:30.709,0.27825,0.25699,0.589582) ",
"ellipse(0:42:22.0130,+41:32:39.489,0.56298,0.33246,162.06169)",
"ellipse(0:42:54.1605,+41:34:29.235,0.27130,0.23529,28.974927)" ]

ra = []
dec = []
for i in ellipses:
    ra.append(i[8:20])
    dec.append(i[21:34])
      
src_pos = []
for i in np.arange(len(ra)):
    src_pos.append( SkyCoord(ra[i], dec[i], unit=(u.hourangle, u.deg)) )

ra_arcsec = []
dec_arcsec = []
for i in src_pos:
    ra_arcsec.append( i.ra.arcsec )
    dec_arcsec.append( i.dec.arcsec )
    
#for i in np.arange(546666-164):
#    ra_arcsec.append( random.uniform(min(ra_arcsec),max(ra_arcsec)) )
#    dec_arcsec.append( random.uniform(min(dec_arcsec),max(dec_arcsec)) )
    
   
#plt.scatter(ra_arcsec,dec_arcsec,s=2,color='blue')
#
#plt.xlabel('ra (arcsecs)')
#plt.ylabel('dec (arcsecs)')
#plt.title('Chandra X-Ray sources')


#print('number of sources:',len(ra_arcsec))


#proper motions and positions of lenses:

v_p = []
d_lens = []
delta_t = []
n_lines = 0
n_crossing = 0
prob_lensing = 0
mags = []

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return [array[idx],idx]

einstein_angle = []
speed = []

with open('outputs1.txt', 'a+') as file:

    # Write the header if the file is empty
    if file.tell() == 0:  # Check if the file is empty
        file.write('lens,u0,A,delta_t,prob_lensing\n')

    for velocity in v: #Rsun/day
        idx = v.index(velocity)

        #for mw lenses:
        phi_v = random.uniform(0,math.pi) #radians
        theta_v = random.uniform(0,math.pi) #for velocities not pointing towards us
        d_kpc = random.choices(d_mw,weights=prob_mw,k=1) #sampling one lens distance , output=[distance]
        
        #for m31 lenses:
        #phi_v = random.uniform(math.pi,2*math.pi) #radians
        #theta_v = random.uniform(math.pi,2*math.pi) #for velocities not pointing towards us
        #d_kpc = random.choices(d_m31,weights=prob_m31,k=1) #sampling one lens distance , output=[distance]
                                                          
        d_kpc = d_kpc[0]       
        d_LS = abs(m31 - d_kpc)                           
                                                          
        v_x = velocity*np.sin(theta_v)*np.cos(phi_v)      
        v_y = velocity*np.sin(theta_v)*np.sin(phi_v)      
        v_t = ( v_x*v_x + v_y*v_y )**0.5 #transverse velocity
        
        v_p = v_t / (4.74*1000*d_kpc)
        v_p_x = (v_t / (4.74*1000*d_kpc))*np.cos(theta_v) 
        
        #print('v_p_x:',v_p_x)
        #speed.append(v_p_x)
        
        #variation of lens positions with time: 
        y = []                                            
        x0 = random.uniform(37500,39500)
        n_t = 300
        dt = np.linspace(start=0,stop=16,num=n_t) #yrs
        lens_ra = x0 + v_p_x*dt #arcsec
        
        m = v_y / v_x

        #print('v_y:',v_y,'v_x:',v_x,'m:', m)
        
        y_c = 148500
        x_c = 38500
        c = y_c - m*x_c
        #print('c:',c)
        #c= 50000
        lens_dec = m*lens_ra  + c  
        # Clip the lens_dec values to the desired range (147500 to 150000)
        lens_dec = np.clip(lens_dec, min(dec_arcsec), max(dec_arcsec))
        n_lines = n_lines + 1      
        
        
        #find distance to closest source:
        closest_src_ra = find_nearest(ra_arcsec,x0)[0]
        closest_src_dec = find_nearest(dec_arcsec,lens_dec[0])[0]
        src_idx_ra = find_nearest(ra_arcsec,x0)[1]
        src_idx_dec = find_nearest(dec_arcsec,lens_dec[0])[1]
        
        initial_diff_dec = abs(closest_src_dec - lens_dec[0])
        initial_diff_ra = abs(closest_src_ra - lens_ra[0])
        
        theta_E = ( mass[idx]/(10**(11.09)) )**0.5 * (d_kpc*m31*10**(-6)/d_LS)**(-0.5)#arcsec
        delta_t = theta_E/v_p #years
        
        #einstein_angle.append(theta_E)
           
        #plt.scatter(lens_ra[-1],lens_dec[-1],s=2,color='red')
        
        if (lens_ra[0] < closest_src_ra < lens_ra[-1] or lens_ra[-1] < closest_src_ra < lens_ra[0]) and (lens_dec[0] < closest_src_dec < lens_dec[-1] or lens_dec[-1] < closest_src_dec < lens_dec[0]) :  #then proceed 
            
            #print('yay!')
            #print('initial ra:',lens_ra[0],'final ra:',lens_ra[-1])
            #print('closest_src_ra:',closest_src_ra)
            #theta_E_rad = theta_E * math.pi/(180 * 3600) #radians
            #r_E = theta_E_rad * m31 
            #x = (d_kpc)/m31 
            delta_t = theta_E/v_p #years
            
            for lens in np.arange(len(lens_dec)): #for every point in lens trajectory
                y_diff = abs(lens_dec[lens]-dec_arcsec[src_idx_dec])
                x_diff = abs(lens_ra[lens]-ra_arcsec[src_idx_ra])
                u0 = np.sqrt(y_diff*y_diff + x_diff*x_diff) / theta_E

                if y_diff < theta_E and x_diff < theta_E:
                    
                    n_crossing = n_crossing + 1
                    
                    #finding magnification:
                    t0 = 0 #time at which lensing occurs 
                    mags = []
                    t_point = dt[lens] #time of lens at the point of crossing, yrs
                    
                    u_t = ( u0**2 + ((t_point-t0)/delta_t)**2 )**(0.5)                          
                    A = (u_t*u_t +2)/(u_t * (u_t*u_t+4)**(0.5))                  
                        
                    prob_lensing = n_crossing / (n_lines*n_t)     
                    file.write(f'\n{n_lines},{u0},{A},{delta_t},{prob_lensing}')
                     
            print(n_lines,prob_lensing,theta_E)

#plt.savefig('/home/sln/summer_2023/plots/srcpositions.png',dpi=300,bbox_inches='tight')
#plt.close()

#print('min and max angle:',min(einstein_angle),max(einstein_angle))
#print('speed:',min(speed),max(speed))
