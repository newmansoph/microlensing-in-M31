# microlensing-towards-M31
A repository containing the code and data required to recreate the results in the paper 'Resolving accretion flows using lensing by white dwarfs in M31' by Sophie Newman and Matthew Middleton (2024). 

For full details of the code and data used, refer to the explanations within the paper.

Within the code folder:
- sim_wd.py simulates the paths of white dwarfs drawn from a 2D mass-velocity distribution as described within the paper and outputs the magnifications, crossing times and impact parameters for the simulated microlensing events 
- paczynski.py uses the outputs of sim_wd.py to create the Paczynski curves for each event
- accretion_disk.py is used to create the flux profile for different bands (XMM-Newton, JWST NIRCam, Swift UVOT, ZTF-g,r,i)
- extended_src.py uses the radii of max flux from accretion_disk.py and the whole flux profile along with the Witt & Mao (1994) equations for magnification to create the magnified flux profiles for each band
- to run these, be in the main directory and type "python ./code/sim_wd.py" for example

Data:
- xamin.csv contains information on the Swift observation towards M31 between September 2006 to March 2023
- m31_catalogue_within_fov.txt contains the positions, summed fluxes, net counts and total exposure times for the Vulic et al. (2016) Chandra observations towards M31 in the mean Swift FOV

Plots:
- Where plots created using the scripts in ./code are saved

Sim results:
- All of the results of sim_wd.py are recorded here
