# Program to measure distance, velocity, and mass of galactic particle

# Import important packages and 'Read' function to read text files
import numpy as np
import astropy.units as u
from ReadFile import Read

# Define function to read text file and determine important particle info
# Paramters are file name, particle type (1=Dark Matter,2=Disk,3=Halo)
def ParticleInfo(filename,partType,partNum):
    
    # Call 'Read' function from ReadFile.py
    time, numParticle, data = Read(filename)

    # Get index for all particles of given type
    index = np.where(data['type']==partType)

    # Read in position components of particle
    x = data[index]['x'][partNum]*u.kpc
    y = data[index]['y'][partNum]*u.kpc
    z = data[index]['z'][partNum]*u.kpc
    
    # Read in velocity components of particle
    vx = data[index]['vx'][partNum]*(u.km/u.s)
    vy = data[index]['vy'][partNum]*(u.km/u.s)
    vz = data[index]['vz'][partNum]*(u.km/u.s)
    
    # Read in mass of each particle 
    # Multiply by 10^10 to convert mass from text file
    m = (data[index]['m'][partNum]*(10**10))*u.solMass
    
    # Calculate magnitude of particle's 3D distance
    distance = np.around(np.sqrt((x**2)+(y**2)+(z**2)),3)
    
    # Calculate magnitude of particle's 3D velocity
    velocity = np.around(np.sqrt((vx**2)+(vy**2)+(vz**2)),3)
    
    # Save particle mass
    mass = m
    
    print(distance,velocity,mass)
    return distance, velocity, mass
ParticleInfo('C:/Users/gpigs_000/400B_Lilly/Homeworks/Homework2/MW_000.txt',2,99)
    
    