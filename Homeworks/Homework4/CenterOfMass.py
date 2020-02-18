# Homework 4
# Center of Mass Position and Velocity
# Jimmy Lilly (2/13/20)

# Import modules
import numpy as np
import astropy.units as u
from ReadFile import Read

# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot
class CenterOfMass:
    
    # Initialize the instance of this Class with the following properties:
    def __init__(self, filename, ptype):
    
        # Read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        # Create array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)
        
        # Store the positions, velocities, and mass of only the particles of the given type
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]
        self.m = self.data['m'][self.index]

    # Function to compute the center of mass position or velocity generically
    def COMdefine(self,a,b,c,m):
    # input: 
    #    array (a,b,c) of positions or velocities and the mass
    # returns: 
    #    3 floats  (the center of mass coordinates)

        # Compute generic COM
        # xcomponent Center of mass
        xProduct = np.multiply(a,m) # x-component times mass of this component
        Acom = np.sum(xProduct)/np.sum(m)
        
        # ycomponent Center of mass
        yProduct = np.multiply(b,m) # y-component times mass of this component
        Bcom = np.sum(yProduct)/np.sum(m)
        
        # zcomponent Center of mass
        zProduct = np.multiply(c,m) # z-component times mass of this component
        Ccom = np.sum(zProduct)/np.sum(m)
        
        return Acom, Bcom, Ccom
    
    # Function to specifically return the center of mass position
    def COM_P(self, delta):                                         
    # Input:                                                                                                           
    #     particle type (1,2,3)                                                                                     
    #     delta (tolerance)                                                                                         
    # Returns: 
    #     One vector, with rows indicating:                                                                                                                                                                            
    #      3D coordinates of the center of mass position (kpc)                                                             

        # Finding center of mass position                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z,self.m)
        
        # Compute the magnitude of the COM position vector.
        RCOM = np.sqrt(XCOM**2+YCOM**2+ZCOM**2)

        # Iterative process to determine the center of mass                                                            

        # Change reference frame to COM frame                                                                          
        # Compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        xNew = self.x - XCOM
        yNew = self.y - YCOM
        zNew = self.z - ZCOM
        RNEW = np.sqrt(xNew**2+yNew**2+zNew**2)

        # Find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        RMAX = (max(RNEW)/2.0)
        
        # Pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        CHANGE = 1000.0

        # Start iterative process to determine center of mass position                                                 
        # Delta is the tolerance for the difference in the old COM and the new one.    
        while (CHANGE > delta):
            
            # Select all particles within the reduced radius (starting from original x,y,z,m)
            index2 = np.where(RNEW <= RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # Compute the center of mass position using                                                                
            # the particles in the reduced radius
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2, y2, z2, m2)
            
            # Compute the new 3D COM position
            RCOM2 = np.sqrt(XCOM2**2+YCOM2**2+ZCOM2**2)

            # Determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)
            # uncomment the following line if you wnat to check this                                                                                               
            #print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            RMAX = RMAX/2.0
            # check this.                                                                                              
            #print ("RMAX", RMAX)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            xNew = self.x - XCOM2
            yNew = self.y - YCOM2
            zNew = self.z - ZCOM2
            RNEW = np.sqrt(xNew**2+yNew**2+zNew**2)

            # Set the center of mass positions to the refined values                                                   
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2 

            # Create a vector to store the COM position with correct units                                                                                                                                                   
            COMP = [XCOM, YCOM, ZCOM]
            
        return np.round(COMP*u.kpc,2)
    
    # Function calculate the center of mass velocity
    def COM_V(self, COMX,COMY,COMZ):
        # Input: 
        #    X, Y, Z positions of the COM
        # Returns: 
        #    3D Vector of COM Velocities
        
        # The max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0*u.kpc

        # Determine the position of all particles relative to the center of mass position
        xV = self.x*u.kpc - COMX
        yV = self.y*u.kpc - COMY
        zV = self.z*u.kpc - COMZ
        RV = np.sqrt(xV**2+yV**2+zV**2)
        
        # Determine the index for those particles within the max radius
        indexV = np.where(RV < RVMAX)

        # Determine the velocity and mass of those particles within the mass radius
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew = self.m[indexV]
        
        # Compute the center of mass velocity using those particles
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew,vynew,vznew,mnew)

        # Create a vector to store the COM velocity
        COMV = [VXCOM, VYCOM, VZCOM]

        # Return the COM vector with appropriate units                                                                                    
        return np.round(COMV*u.km/u.s,2)
    
# ANSWERING QUESTIONS
#######################

# Create a center of mass object for the MW, M31 and M33
MWCOM = CenterOfMass('C:/Users/Jimmy/400B_Lilly/MW_000.txt', 2)
M31COM = CenterOfMass('C:/Users/Jimmy/400B_Lilly/M31_000.txt', 2)
M33COM = CenterOfMass('C:/Users/Jimmy/400B_Lilly/M33_000.txt', 2)

# MW: store the COM position and velocity
MW_COMP = MWCOM.COM_P(0.1)
MW_COMV = MWCOM.COM_V(MW_COMP[0],MW_COMP[1],MW_COMP[2])
print("MW COM Position " + str(MW_COMP))
print("MW COM Velocity " + str(MW_COMV))
print('\n')

# M31: store the COM position and velocity
M31_COMP = M31COM.COM_P(0.1)
M31_COMV = M31COM.COM_V(M31_COMP[0],M31_COMP[1],M31_COMP[2])
print("M31 COM Position " + str(M31_COMP))
print("M31 COM Position " + str(M31_COMV))
print('\n')

# M33: store the COM position and velocity
M33_COMP = M33COM.COM_P(0.1)
M33_COMV = M33COM.COM_V(M33_COMP[0],M33_COMP[1],M33_COMP[2])
print("M33 COM Position " + str(M33_COMP))
print("M33 COM Velocity" + str(M33_COMV))
print('\n')

# Determine relative separation between MW and M31
MW_M31_x_diff = MW_COMP[0]-M31_COMP[0]
MW_M31_y_diff = MW_COMP[1]-M31_COMP[1]
MW_M31_z_diff = MW_COMP[2]-M31_COMP[2]
MW_M31_separation = np.round(np.sqrt(MW_M31_x_diff**2+MW_M31_y_diff**2+MW_M31_z_diff**2),2)
print("Relative separation between MW and M31 = " + str(MW_M31_separation))

# Determine magnitude of relative velocity between MW and M31
MW_M31_vx_diff = MW_COMV[0]-M31_COMV[0]
MW_M31_vy_diff = MW_COMV[1]-M31_COMV[1]
MW_M31_vz_diff = MW_COMV[2]-M31_COMV[2]
MW_M31_velocity = np.round(np.sqrt(MW_M31_vx_diff**2+MW_M31_vy_diff**2+MW_M31_vz_diff**2),2)
print("Relative velocity between MW and M31 = " + str(MW_M31_velocity))
print('\n')

# Determine relative separation between M31 and M33
M31_M33_x_diff = M31_COMP[0]-M33_COMP[0]
M31_M33_y_diff = M31_COMP[1]-M33_COMP[1]
M31_M33_z_diff = M31_COMP[2]-M33_COMP[2]
M31_M33_separation = np.round(np.sqrt(M31_M33_x_diff**2+M31_M33_y_diff**2+M31_M33_z_diff**2),2)
print("Relative separation between M31 and M33 = " + str(M31_M33_separation))

# Determine magnitude of relative velocity between M31 and M33
M31_M33_vx_diff = M31_COMV[0]-M33_COMV[0]
M31_M33_vy_diff = M31_COMV[1]-M33_COMV[1]
M31_M33_vz_diff = M31_COMV[2]-M33_COMV[2]
M31_M33_velocity = np.round(np.sqrt(M31_M33_vx_diff**2+M31_M33_vy_diff**2+M31_M33_vz_diff**2),2)
print("Relative velocity between M31 and M33 = " + str(M31_M33_velocity))

"""
# Question 4
#######################
 Q: Given that M31 and the MW are about to merge, why is the iterative process 
    to determine the COM important?
 A: The iterative process is important because, as the galaxies continue to get 
    closer, the dynamics of their merger will be determined by the common center 
    of mass between the galaxies rather than their individual centers of mass 
    which we have defined thus far in this program.
#######################
"""
