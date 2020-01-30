# Program to read in text file and assign AstroPy units to data in the file

# Import important packages (NumPy and AstroPy)
import numpy as np
import astropy.units as u

# Define function to read file and extract data
def Read(filename):
    
    # Open text file
    file = open(filename,'r')
    
    # Read 1st row of text file
    line1 = file.readline()
    # Split data in 1st row (delimiter is space)
    label1, value1 = line1.split() 
    # Save time data from 1st row in 10^6 yrs
    time = float(value1)*u.Myr
    
    # Read 2nd row
    line2 = file.readline()
    # Split data in 2nd row (delimiter is space)
    label2, value2 = line2.split()
    # Save data from 2nd row as number of particles
    numParticles = float(value2)
    
    # Close text file that was opened
    file.close() 
    
    # Store remaining data from text file in an array
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
    # Return time of snapshot, # of particles, and array of other data
    return time, numParticles, data