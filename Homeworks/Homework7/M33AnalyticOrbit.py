# Homework 7
# Analytical Orbit of M33
# Jimmy Lilly (4/3/20)

# Import modules
import numpy as np
import astropy.units as u
from ReadFile import Read
import astropy.constants as const
from CenterOfMass2 import CenterOfMass
from GalaxyMass import ComponentMass

# Class to determine the acceleration M33 feels from M31 and integrate its 
# current position and velocity forwards in time
class M33AnalyticalOrbit:
    
    # Initialize the instance of this Class with the following properties:
    def __init__(self, filename):
    
        # Define gravitational constant as global variable
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value

        # Define centers of mass for M31 and M33 by calling CoM function
        M31COM = CenterOfMass('C:/Users/Jimmy/400B_Lilly/M31_000.txt', 2)
        M33COM = CenterOfMass('C:/Users/Jimmy/400B_Lilly/M33_000.txt', 2)
        
        M31_COMP = M31COM.COM_P(0.1,2) # delta = 0.1, VolDec = 2
        M31_COMV = M31COM.COM_V(M31_COMP[0],M31_COMP[1],M31_COMP[2])
        
        M33_COMP = M33COM.COM_P(0.1,4) # delta = 0.1, VolDec = 4
        M33_COMV = M33COM.COM_V(M33_COMP[0],M33_COMP[1],M33_COMP[2])
        
        # Calculate relative position and velocity of M31 and M33
        self.r0 = (M33_COMP-M31_COMP).value
        self.v0 = (M33_COMV-M31_COMV).value
        
        # Assign filename to variable
        self.filename = filename
        
        # Scale lengths and masses for each component of M31
        self.rdisk = 5.0 # in kpc
        self.Mdisk = 1e12*ComponentMass('C:/Users/Jimmy/400B_Lilly/M31_000.txt',2).value # in Msun
        
        self.rbulge = 1.0 # in kpc
        self.Mbulge = 1e12*ComponentMass('C:/Users/Jimmy/400B_Lilly/M31_000.txt',3).value # in Msun
        
        self.rhalo = 61.5 # in kpc
        self.Mhalo = 1e12*ComponentMass('C:/Users/Jimmy/400B_Lilly/M31_000.txt',1).value # in Msun
        
    # Function to calculate accel. of halo and bulge using Hernquist profile
    def HernquistAccel(self, M, r_a, r):
        # Inputs:
        #    M = mass of galaxy component (in Msun)
        #    r_a = scale length of galaxy component (in kpc)
        #    r = relative position vector of galaxy
        # Returns:
        #    Acceleration vector from a Hernquist potential
        
        # Magnitude of the relative position vector
        rmag = np.sqrt(r[0]**2+r[1]**2+r[2]**2)
        
        # Prefactor 
        prefactor = -self.G*M/rmag*((r_a+rmag)**2)
        
        # Calculate acceleration of component
        a = prefactor*self.r0
        
        # Return acceleration vector of bulge or halo
        return a
    
    # Function to calculate accel. of disk using Miyamoto-Nagai 1975 profile
    def MiyamotoNagaiAccel(self, M, rd, r):
        # Inputs:
        #    M = mass of galaxy's disk (in Msun)
        #    rd = disk radius (in kpc)
        #    r = relative position vector of galaxy
        # Returns:
        #    Acceleration vector of galaxy disk from Miyamoto-Nagai profile
        
        # Define z scale height for disk
        zd = self.rdisk/5.0
        
        # Define R and B parameters from Miyamoto-Nagai 1975 profile
        R = np.sqrt(r[0]**2+r[1]**2)
        B = rd + np.sqrt(r[2]**2+zd**2)
        
        # Factors to multiple each acceleration component by
        multiplier = np.array([1,1,B/np.sqrt(r[2]**2+zd**2)])
        
        # Define common prefactor for each acceleration component
        prefactor = (-self.G*M/((R**2+B**2)**1.5))*multiplier
        
        # Calculate components of acceleration vector
        a = np.multiply(np.array(r),prefactor)
        
        # Return acceleration vector of disk
        return a
    
    # Function to sum acceleration vectors of each galaxy component
    def M31Accel(self,r):
        #    r = relative position vector of galaxy
        # Returns:
        #    Acceleration vector of galaxy disk from Miyamoto-Nagai profile
        
        # Calculate acceleration vector of each component
        # Use functions defined above
        a_halo = self.HernquistAccel(self.Mhalo,self.rhalo,self.r0)
        a_bulge = self.HernquistAccel(self.Mbulge,self.rbulge,self.r0)
        a_disk = self.MiyamotoNagaiAccel(self.Mdisk,self.rdisk,self.r0)
        
        # Add halo and bulge accelerations (np.add can only add 2 lists)
        sum1 = np.add(a_halo,a_bulge)
        
        # Add halo+bulge to disk
        acceleration = np.add(sum1,a_disk)
        
        # Return total acceleration vector of M31
        return acceleration
    
    # Function to integrate the acceleration vector 
    # (forward: dt > 0 or backward: dt < 0)
    def LeapFrog(self, dt, r, v):
        # Inputs:
        #    dt: time interval for integration
        #    r: relative position vector of M33 and M31
        #    v: relative velocity vector of M33 and M31
        # Returns:
        #    Updated position and velocity vectors
        
        # Predict the position at the next half timestep
        rhalf = r + v*dt/2
        
        # Compute the velocity at the next timestep
        vnew = v + self.M31Accel(rhalf)*dt
        
        # Compute the position at the next timestep
        rnew = rhalf + vnew*dt*.5
        
        # Return the new position and velcoity vectors
        return rnew, vnew
    
    # Function to loop over LeapFrog to solve the EoM and compute
    # the future orbit of M33 for 10 Gyr into the future
    def OrbitIntegration(self, t0, dt, tmax):
        # Inputs:
        #    t0 = time to start assessing orbit
        #    dt = time interval
        #    tmax = time to stop assessing orbit
        # Returns:
        #    Future orbit of M31

        # Initialize the time to the input starting time
        t = t0
        
        # Initialize an empty array of size:  rows-int(tmax/dt)+2, columns-7
        orbit = np.zeros((int(tmax/dt)+2,7))
        
        # Initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        
        # This above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        # Initialize a counter for the orbit
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # Start the integration (advancing in time steps and computing LeapFrog at each step)
        while (t<tmax):  # as long as t has not exceeded the maximal time 
            
            # Advance the time by one timestep, dt
            t += dt
           
            # Store the new time in the first element of the ith row
            orbit[i][0] = t
            
            # Advance the position and velocity using the LeapFrog scheme
            rnew, vnew = self.LeapFrog(dt, orbit[i-1,1:4], orbit[i-1,4:7])
            print(rnew)
            print(vnew)
            print("Loop "+str(i)+" has completed")
            print('\n')
            
            # Store new position vector into the columns with indexes 1,2,3
            # of the ith row of orbit
            orbit[i, 1:4] = rnew 
            
            # Store new velocity vector into columns with indexes 4,5,6
            # of the ith row of orbit
            orbit[i, 4:7] = vnew
            
            # Update counter i , where i is keeping track of the number of rows 
            # (i.e. the number of time steps)
            i += 1
        
        # Write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))

"""
# Define name/destination of text file
file = 'C:/Users/Jimmy/400B_Lilly/Homeworks/Homework7/OrbitIntegration.txt'

# Initialize instance of class
Orbit = M33AnalyticalOrbit(file)

# Call OrbitIntegration function
Orbit.OrbitIntegration(0.0,.01,10.0)
"""