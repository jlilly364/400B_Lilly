# Homework 6
# Plotting orbits of COM's of MW, M31, and M33
# Jimmy Lilly (3/6/20)
# Worked with Sean Cunningham, Madison Walder, Ryan Webster

# Import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# Import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# Import modules from previous homeworks
from ReadFile import Read

# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass2 import CenterOfMass

# Function to loop over snapshots to find COM pos and vel. vs. time
def OrbitCOM(galaxy,start,end,n):
    # Inputs:
    #     galaxy = name of galaxy (MW, M31, or M33)
    #     start = number of first snapshot
    #     end = number of last snapshot
    #     n = intervals at which to return COM
          
    # Returns: 
    #     COM_P = position of center of mass of galaxy
    #     COM_V = velocity of center of mass of galaxy
    
    # compose the filename for output
    fileout = "C:/Users/Jimmy/400B_Lilly/Homeworks/Research Assignments/3/Orbit_{0}.txt".format(galaxy)
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    if galaxy == 'M33':
        delta = 0.1
        VolDec = 4
    else:
        delta = 0.1
        VolDec = 2
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids = np.arange(start,end,n)
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([len(snap_ids), 7])
    
    # Loop to calculate COM pos. and vel. at each snapshot
    for i, snap_id in enumerate(snap_ids):
        
        # compose the data filename (be careful about the folder)
        # Add string of filenumber to value 000
        ilbl = '000' + str(snap_id)
        
        # Remove all but last 3 digits of string
        ilbl = ilbl[-3:]
        
        # Assign filename based on snapshot and galaxy inputs
        filename = 'C:/Users/Jimmy/Downloads/' + '%s/'%(galaxy) + '%s_'%(galaxy) + ilbl + '.txt'
        
        # Initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(filename, 2)

        # Store the COM pos and vel. Remember that now COM_P required VolDec
        COM_P = COM.COM_P(delta,VolDec)
        COM_V = COM.COM_V(COM_P[0],COM_P[1],COM_P[2])
    
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store: a[i] = var1, *tuple(array1)
        orbit[i] = COM.time.value/1000, *tuple(COM_P.value), *tuple(COM_V.value)
        
        # print snap_id to see the progress
        print('Calculating COM info for '+str(snap_id))

    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
    
# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 
# Uncomment to make new OrbitCOM text files

OrbitCOM("MW",0,801,1)

