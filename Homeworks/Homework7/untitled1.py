"""
Creation date: Thu Apr  2 19:20:33 2020
Author: Jimmy Lilly (www.github.com/jlilly364)

Program Objective: 
"""
import astropy.units as u
import astropy.constants as const
import numpy as np
from GalaxyMass import ComponentMass

G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value

r = [1,2,3]
rhalo = 61.5 # in kpc
Mhalo = 1e12*ComponentMass('C:/Users/Jimmy/400B_Lilly/M31_000.txt',1).value # in Msun

rdisk = 5.0
Mdisk = 1e12*ComponentMass('C:/Users/Jimmy/400B_Lilly/M31_000.txt',2).value # in Msun
print("Mdisk = "+str(Mdisk))
zd = rdisk/5.0

# Function to calculate accel. of halo and bulge using Hernquist profile
def HernquistAccel(M, r_a, r):
    # Inputs:
    #    M = mass of galaxy component (in Msun)
    #    r_a = scale length of galaxy component (in kpc)
    #    r = relative position vector of galaxy
    # Returns:
    #    Acceleration vector from a Hernquist potential
    
    # Magnitude of the relative position vector
    rmag = np.sqrt(r[0]**2+r[1]**2+r[2]**2)
    print("rmag = "+str(rmag))
    
    # Prefactor 
    A = -G*M
    B = rmag*((r_a+rmag)**2)
    prefactor = A/B
    print("prefactor = "+str(prefactor))
    
    # Calculate acceleration of component
    a = [x*prefactor for x in r]
    print(a)
    
    # Return acceleration vector of bulge or halo
    return a

# Function to calculate accel. of disk using Miyamoto-Nagai 1975 profile
def MiyamotoNagaiAccel(M, rd, r):
    # Inputs:
    #    M = mass of galaxy's disk (in Msun)
    #    rd = disk radius (in kpc)
    #    r = relative position vector of galaxy
    # Returns:
    #    Acceleration vector of galaxy disk from Miyamoto-Nagai profile
    
    # Define z scale height for disk
    zd = rdisk/5.0
    print("zd = "+str(zd))
    
    # Define R and B parameters from Miyamoto-Nagai 1975 profile
    R = np.sqrt(r[0]**2+r[1]**2)
    print("R = "+str(R))
    B = rd + np.sqrt(r[2]**2+zd**2)
    print("B = "+str(B))
    
    # Factors to multiple each acceleration component by
    multiplier = np.array([1,1,B/np.sqrt(r[2]**2+zd**2)])
    print(multiplier)
    
    # Define common prefactor for each acceleration component
    prefactor = (-G*M/((R**2+B**2)**1.5))*multiplier
    print(prefactor)
    
    # Calculate components of acceleration vector
    a = np.multiply(np.array(r),prefactor)
    print(a)
    
    # Return acceleration vector of disk
    return a

#HernquistAccel(Mhalo,rhalo,r)
MiyamotoNagaiAccel(Mdisk,rdisk,r)