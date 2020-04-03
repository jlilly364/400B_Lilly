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
    fileout = "C:/Users/Jimmy/400B_Lilly/Homeworks/Homework6/Orbit_{0}.txt".format(galaxy)
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    if galaxy == 'M33':
        delta = 0.1
        VolDec = 4
    else:
        delta = 0.1
        VolDec = 5
    
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
    """
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
    """

# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 
# Uncomment to make new OrbitCOM text files
"""
OrbitCOM("MW",0,800,5)
OrbitCOM("M31",0,800,5)
OrbitCOM("M33",0,800,5)
"""

# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
MW_file = "C:/Users/Jimmy/400B_Lilly/Homeworks/Homework6/Orbit_MW.txt"
data_MW = np.genfromtxt(MW_file,dtype=None,names=True)

M31_file = "C:/Users/Jimmy/400B_Lilly/Homeworks/Homework6/Orbit_M31.txt"
data_M31 = np.genfromtxt(M31_file,dtype=None,names=True)

M33_file = "C:/Users/Jimmy/400B_Lilly/Homeworks/Homework6/Orbit_M33.txt"
data_M33 = np.genfromtxt(M33_file,dtype=None,names=True)

# Function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit
def VectorDiff(data1,data2):
    # Inputs:
    #     vector1: first vector in pair
    #     vector2: second vector in pair
    # Returns:
    #    diffmag_P = magnitude of difference between position vectors
    #    diffmag_V = magnitude of difference between velocity vectors

    # Calculate differences in position coordinates
    xdiff = np.subtract(data1['x'],data2['x'])
    ydiff = np.subtract(data1['y'],data2['y'])
    zdiff = np.subtract(data1['z'],data2['z'])
    
    # Calculate differences in velocity coordinates
    vxdiff = np.subtract(data1['vx'],data2['vx'])
    vydiff = np.subtract(data1['vy'],data2['vy'])
    vzdiff = np.subtract(data1['vz'],data2['vz'])

    # Calculate magnitudes of vector differences
    diffmag_P = np.sqrt(xdiff**2+ydiff**2+zdiff**2)
    diffmag_V = np.sqrt(vxdiff**2+vydiff**2+vzdiff**2)

    # Save time from each line of data
    time = data1['t']

    return time, diffmag_P, diffmag_V

# Determine the magnitude of the relative position and velocities 
    
# of MW and M31
timeMW_M31, MW_M31_diffP, MW_M31_diffV = VectorDiff(data_MW,data_M31)

# of M31 and M33
timeM31_M33, M31_M33_diffP, M31_M33_diffV = VectorDiff(data_M31,data_M33)

# Function to plot orbits
def OrbitPlot(time1,position1,time2,position2):
    # Plot the Orbit of the galaxies 
    #################################
    fig,ax = plt.subplots(figsize=(10,8))
    
    # Plot separations vs. time
    #plt.plot(timeMW_M31,MW_M31_diffP,linewidth=5,label='MW-M31')
    #plt.semilogy(timeMW_M31,MW_M31_diffP,linewidth=5,label='MW-M31')
    plt.plot(time1,position1,linewidth=5,label='M31-M33 Simulation',linestyle='dashed')
    plt.plot(time2,position2,linewidth=5,label='M31-M33 Analytic')
    
    # Give labels and title
    plt.xlabel("Time (Gyr)",fontsize=22)
    plt.ylabel("Separation (kpc)",fontsize=22)
    plt.title("Relative Position between M31-M33",fontsize=24)
    
    # Plot legend and show plot
    legend = ax.legend(fontsize='x-large')
    
    # Save plot
    sep_file = 'C:/Users/Jimmy/400B_Lilly/Homeworks/Homework7/Separation.png'
    
    # Uncomment below to save file
    plt.savefig(sep_file)

# Function to plot relative velocities
def VelocityPlot(time1,velocity1,time2,velocity2):
    # Plot the orbital velocities of the galaxies 
    #################################
    fig,ax = plt.subplots(figsize=(10,8))
    
    # Plot separations vs. time
    plt.plot(time1,velocity1,linewidth=5,label='M31-M33 Simulation',linestyle='dashed')
    plt.plot(time2,velocity2,linewidth=5,label='M31-M33 Analytic')
    
    # Give labels and title
    plt.xlabel("Time (Gyr)",fontsize=22)
    plt.ylabel("Relative Velocity (km/s)",fontsize=22)
    plt.title("Relative Velocity between M31-M33",fontsize=24)
    
    # Plot legend and show plot
    legend = ax.legend(fontsize='x-large')
    
    # Save plot
    #vel_file = 'C:/Users/Jimmy/400B_Lilly/Homeworks/Homework6/LG_velocity.png'
    # Uncomment below to save file
    #plt.savefig(vel_file)

"""
Answers to questions:

    1. Q: How many close encounters will the MW and M31 experience in the 
          future?
       A: The MW and M31 will have 3 close encounters in the future
          (~3.97,5.87,6.64 billion years)
    
    2. Q: How is the time evolution of the separation and relative velocity
          related?
       A: When the separation is least (gravitational force between them is 
          greatest), the relative velocity is greatest. When the separation is 
          greatest (gravitational force between them is least), the relative 
          velocity is least. They boomerang around each other until merging.
          
    3. Q: When do M31 and the MW merge? What happen's to M33's orbit as they
          merge?
       A: M31 and the MW seem to merge in 6.64 billion years. When they merge, 
          M33 orbits the remnant galaxy with a regular period and average 
          relative velocity of 250 km/s. Over the next ~6 billion years,
          M33 falls into (gets closer to) the remnant galaxy.
"""