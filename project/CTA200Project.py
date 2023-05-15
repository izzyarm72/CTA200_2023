'''
    Simple catalog

    Modified from myexample.py for CTA200 project
'''
#importing necessary packages for section 1
import numpy as np
from speedystar import starsample
from speedystar.eject import Hills
from astropy import units as u
from speedystar.utils.mwpotential import MWPotential 

'''
    Create an ejection catalog of HVSs and propagate it through varying MW potentials 
'''

#------------------------
'''
Section 1: creating, propagating and perfomring photometry on the star sample
'''


'''
    Create ejection catalog
'''

# Initialize an ejection model, i.e. how the spatial and velocity 
#distribution of the stars will be sampled
#using an ejection rate one order of magnitue larger than default 1e-4
ejectionmodel = Hills(rate = 1e-3/u.yr)

# Eject a sample of stars from Sgr A*. 
mysample = starsample(ejectionmodel, name='CTA200 Q2: MW Potential')

# Save ejection sample 
#this sample will be propagated under different conditions 
mysample.save('./MWP_ejection.fits')

'''
    Propagate ejection catalogue through the galaxy.
    Propatate once through default potential and once through a potential 
    with a greater halo mass. 
'''

#list of different MW halo mass values
#0.76 is defualt, testing order of magnitude higher & lower 
M_s = [0.076, 0.76, 1.76]

for i, mass in enumerate(M_s):
    
    #making the sample
    mysample = starsample('./MWP_ejection.fits')
    
    # use mwpotential from https://arxiv.org/abs/2108.01100
    Potential = MWPotential(Ms=mass, rs=24.8, c=1., T=True)
    
    
    #propagate 
    mysample.propagate(potential = Potential)
    
    #save propagated sample 
    mysample.save('./MWP_propagated%i.fits'%i)


    '''
    Calculate the apparent magnitudes of each star 
    '''

    #Load in MWP_propagated samples that was just generated
    mysample = starsample('./MWP_propagated%i.fits'%i)
    
    #Selecting HVS with velocities greater than the escape velocity to infinity 
    #These are the HVS we are intereste in 
    #Will also speed up program run time
    idx = (mysample.GCv>mysample.Vesc) 
    mysample.subsample(np.where(idx)[0])
    
    
    #Magnitudes are excincted by Milky Way dust along the line of sight. 
    #Before first use, dust maps must be downloaded, see docstring of
    #speedystar.fetch_dust()

    #change this line before submitting so that it can be run by others
    #uncomment 
    mysample.fetch_dust(path='/path/to/wherever/') 
    
    mysample.config_dust('./../dust-map-3d.h5')
    mysample.photometry()

    #now save photometry 
    mysample.save('./MWP_photometry%i.fits'%i)

#----------------------

'''
    Section 2: Explore the differences between the three samples.
    Will extract information from the photoemtry fits files, and plot results
'''


#import necessary pacakges for section 2
from astropy.io import fits
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})


'''
    Read in the fits files
'''

def load_fits(file):
    """Loads in a fits file at extracts data 
    -------
    Parameters:
        file: string 
            path to the file
    -------
    Reutrns:
        GC_v: array
            Galactocentric Velocity of the observable HVSs 
        V_esc: array
            Escape velocity to inifinity for given HVS position 
        GC_dist: array 
            Galactocentric total distance
        obs_HVS: array 
            indices of HVSs with GaiaG < 20.7
            
    """
    hdul = fits.open(file)
    HVS_pop = np.array(hdul[1].data)
    
    #read in magnitudes  
    GaiaG = HVS_pop['Gaia_G']
    
    #selecting HVS that are detectable by Gaia  
    obs_HVS = np.where((GaiaG < 20.7))[0]
    
    #extract desired infromation    
    V_esc = HVS_pop['Vesc'][obs_HVS]
    GC_v = HVS_pop['GCv'][obs_HVS]
    GC_dist = HVS_pop['GCdist'][obs_HVS]
    
    #close the fits file 
    hdul.close()
    
    #return the output 
    return GC_v, V_esc, GC_dist, obs_HVS


'''
    Extract the useful infromation from the .fits files 
    Plots 'Vesc' vs 'GCv'
    Plot 'HVS number' vs 'Gcdist' to understand how far HVS are 
'''

#loading in the three samples 
HVS0 = load_fits('./MWP_photometry0.fits')
HVS1 = load_fits('./MWP_photometry1.fits')
HVS2 = load_fits('./MWP_photometry2.fits')


#now plot  
names = [HVS0, HVS1, HVS2]

#plotting Vesc vs GCv
fig1 = plt.figure(figsize = (8,6))
for i, mass in enumerate(M_s):    
    #plot the three cases 
    plt.plot(names[i][0], names[i][1], '*', label='$M_s$ = %.3f'%mass)
#add labels 
plt.xlabel('Galactocentric Velocity (km/s)', fontsize=15)
plt.ylabel('Escape Velocity (km /s)', fontsize=15)
plt.title('Observable HVS for different Halo Masses')
plt.grid()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#saving the plot
plt.savefig('VescVsGCv.png',  bbox_inches='tight') 
plt.savefig('VescVsGCv.pdf', bbox_inches='tight') 
plt.show()
  

#plotting distances of observable HVS 
fig2 = plt.figure(figsize = (8,6))
for i, mass in enumerate(M_s): 
    #plot the three cases
    plt.plot(names[i][3], names[i][2], '*', label='$M_s$ = %.3f'%mass)
#adding labels 
plt.xlabel('HVS number', fontsize=15)
plt.ylabel('Galactocentric Distance (kpc)', fontsize=15)
plt.title('Distance of Observable HVSs')
plt.grid()
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#saving the plot 
plt.savefig('GCdist.png', bbox_inches='tight')
plt.savefig('GCdist.pdf', bbox_inches='tight')
plt.show()

#----------------------


'''
    Future steps are to investigate changing the axis ration of MW halo 
    and changing the mass of the central MBH 
    to see how those parameters influence the HVS population
'''

#------------------










