# Program to measure distance, velocity, and mass of galactic particle

# Import important packages and 'Read' function to read text files
import numpy as np
import astropy.units as u
from ReadFile import Read

# Define function to read text file and determine important particle info
# Paramters are file name, particle type (1=Dark Matter,2=Disk,3=Halo),
# and particle number
def ParticleInfo(filename,partType,partNum):
    
    # Call 'Read' function from ReadFile.py
    time, numParticle, data = Read(filename)

    # Get index for all particles of given type
    index = np.where(data['type']==partType)

    # Read in position components of particle (in kpc)
    x = data[index]['x'][partNum]*u.kpc
    y = data[index]['y'][partNum]*u.kpc
    z = data[index]['z'][partNum]*u.kpc
    
    # Read in velocity components of particle (in km/s)
    vx = data[index]['vx'][partNum]*(u.km/u.s)
    vy = data[index]['vy'][partNum]*(u.km/u.s)
    vz = data[index]['vz'][partNum]*(u.km/u.s)
    
    # Read in mass of particle 
    # Multiply by 10^10 to convert mass from text file
    mass = (data[index]['m'][partNum]*(10**10))*u.solMass
    
    # Calculate magnitude of particle's 3D distance
    distance = np.around(np.sqrt((x**2)+(y**2)+(z**2)),3)
    
    # Convert particle's 3D distance to lightyears
    distance_lyr = distance.to(u.lyr)
    
    # Calculate magnitude of particle's 3D velocity
    velocity = np.around(np.sqrt((vx**2)+(vy**2)+(vz**2)),3)
    
    # Print calculated particle properties
    print("3D Distance of 100th disk particle: {0}".format(distance))
    print("3D Distance of 100th disk particle: {0}".format(distance_lyr))
    print("3D Velocity of 100th disk particle: {0}".format(velocity))
    print("Mass of 100th disk particle: {0}".format(mass))
    
    return distance, distance_lyr, velocity, mass

# Call function to get info for 100th disk particle
ParticleInfo('C:/Users/gpigs_000/400B_Lilly/Homeworks/Homework2/MW_000.txt',2,99)
    
    