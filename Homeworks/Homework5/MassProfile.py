# Homework 5
# Determine galaxy mass distribution and rotation curve
# Jimmy Lilly (2/24/20)
# Worked with: primarily -> Sean Cunningham, Ryan Webster
#              also -> Mackenzie James, Madison Walder

# Import relevant modules
import numpy as np
import astropy.units as u
from ReadFile import Read
from CenterOfMass import CenterOfMass
from astropy.constants import G
import matplotlib.pyplot as plt

# Convert G to more useful units
G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)

# Create MassProfile class to calculate related quantities for each galaxy
class MassProfile:
    
    # Initialize the instance of this Class with the following properties:
    def __init__(self, galaxy, snap):
        # Inputs:
        #    galaxy = name of galaxy
        #    snap = snapshot number
        # Returns:
        #    Important data from snapshot text files
        
        # Add string of filenumber to value 000
        ilbl = '000' + str(snap)
        
        # Remove all but last 3 digits of string
        ilbl = ilbl[-3:]
        
        # Assign filename based on snapshot and galaxy inputs
        self.filename = 'C:/Users/Jimmy/400B_Lilly/' + '%s_'%(galaxy) + ilbl + '.txt'
    
        # Read data in the given file using Read
        self.time, self.total, self.data = Read(self.filename)

        # Save galaxy name as global variable
        self.gname = galaxy       

        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
        self.m = self.data['m']                                                                                 
    
    # Function to compute mass within given radius of COM for specified galaxy and galaxy component
    def MassEnclosed(self,ptype,radii):
        # Inputs:
        #    ptype = type of particle (1=halo, 2=disk, 3=bulge)
        #    radii = array of radii within given radius of COM position
        # Returns:
        #    masses = array of masses to calculate mass profile with later
        
        # Create array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)
        
        # Store the positions, velocities, and mass of only the particles of the given type
        x = self.x[self.index]
        y = self.y[self.index]
        z = self.z[self.index]
        m = self.m[self.index]
        
        # Define empty list of enclosed masses that's same size as radii list
        M_Enc = []
        
        # Determine COM position of galaxy
        COM = CenterOfMass(self.filename,ptype)
        COM_P = COM.COM_P(0.1)
        
        # Define components of COM
        XCOM = COM_P[0]
        YCOM = COM_P[1]
        ZCOM = COM_P[2]
        
        # Determine COM radius
        RCOM = np.sqrt((x-XCOM)**2+(y-YCOM)**2+(z-ZCOM)**2)
        
        # Measure mass enclosed within each radius
        for i in range(0,len(radii)):
            
            # Find indices of particles enclosed within RCOM 
            new_index = np.where(RCOM < radii[i])
            
            # Fill out M_Enc array accordingly
            M_Enc.append(np.sum(m[new_index])*1e10)
            
        # Uncomment to test code
        #print(M_Enc*u.Msun)
        return M_Enc*u.Msun
    
    # Function to compute total enclosed mass for a galaxy (bulge+disk+halo)
    def MassEnclosedTotal(self,radii):
        # Inputs:
        #     radii = 1D list of radii (in kpc)
        # Returns:
        #     array of masses representing total enclosed mass at each radius
        
        # For MW and M31 include bulge mass in total
        if self.gname != 'M33':
            M_Enc_Halo = self.MassEnclosed(1,radii)
            M_Enc_Disk = self.MassEnclosed(2,radii)
            M_Enc_Bulge = self.MassEnclosed(3,radii)
            M_Enc_Total = M_Enc_Bulge + M_Enc_Disk + M_Enc_Halo
            
        # For M33 exclude bulge mass in total since it doesn't have one!
        else:
            M_Enc_Halo = self.MassEnclosed(1,radii)
            M_Enc_Disk = self.MassEnclosed(2,radii)
            M_Enc_Total = M_Enc_Disk + M_Enc_Halo
        
        # Uncomment to test code
        #print(M_Enc_Total)
        return M_Enc_Total
    
    # Function that computes mass enclosed for a Hernquist profile
    def HernquistM(self, radii, a, Mhalo):
        # Inputs:
        #    r: distance from Galactic center (in kpc)
        #    a: scale radius (in kpc)
        #    Mhalo: total dark matter halo mass (in units of 1e12 Msun)
        # Returns:
        #    Dark matter halo mass at given radius
        
        # Define 
        numerator = Mhalo*1e12*(r**2)
        denominator = (a+r)**2
        HernquistMass = np.round(numerator/denominator,2)*u.Msun
        
        # Uncomment to test function
        #print(HernquistMass)
        return HernquistMass
    
    # Function to calculate the circular velocity at each radius (assuming spherical symmetry)
    def CircularVelocity(self, ptype, radii):
        # Inputs:
        #     ptype = type of particle (1=bulge, 2=disk, 3=halo)
        #     radii = 1D list of radii (in kpc)
        # Returns:
        #     array of circular speeds (in km/s)
        
        # Calculate mass enclosed for specific particle type at each radius
        mass = self.MassEnclosed(ptype, radii)

        # Calculate circular velocity of each component at each radius
        velocity = np.sqrt(G*mass/radii)
            
        return velocity

    # Function to calculate total circular velocity at each radius
    def CircularVelocityTotal(self,radii):
        # Inputs: 
        #     radii = 1D list of radii (in kpc)
        # Returns:
        #     array of circular velocities (in km/s) representing the total Vcirc
        #     created by all the galaxy components (bulge+disk+halo) at each radius
        
        # Define total mass within radius using MassEnclosedTotal
        mass_total = self.MassEnclosedTotal(radii)
        
        # Calculate total circular velocity at each radius
        velocity_total = np.sqrt(G*mass_total/radii)
            
        return velocity_total
    
    # Function to compute circular speed with Hernquist mass profile
    def HernquistVCirc(self, radii, a, Mhalo):
        # Inputs:
        #    r: distance from Galactic center (in kpc)
        #    a: scale radius (in kpc)
        #    Mhalo: total dark matter halo mass (in units of 1e12 Msun)
        # Returns:
        #    Dark matter halo mass at given radius
        
        # Calculate Hernquist mass for given radius
        mass = self.HernquistM(radii, a, Mhalo)

        # Calculate Hernquist velocity at each radius
        velocity_Hernq = np.sqrt(G*mass/radii)
            
        return np.round(velocity_Hernq,2)
    
# Generate array of radii from 0.1 to 30.5 kpc
r = np.arange(0.1,30.5,.25)*u.kpc

# Uncomment below to make plots for MW
"""
# Call class for snapshot 0 of MW
MW = MassProfile('MW',0)

# Initialize MW scale radius and halo mass (from HW 3)
MW_a = 61.5*u.kpc
MW_Mhalo = 1.975

# Plot mass profile of Milky Way (each component and total)
plt.figure(1)
plt.semilogy()
plt.plot(r,MW.MassEnclosed(1,r),label='Halo Mass',linestyle='dashed')
plt.plot(r,MW.MassEnclosed(2,r),label='Disk Mass',linestyle=':')
plt.plot(r,MW.MassEnclosed(3,r),label='Bulge Mass',linestyle='dashdot')
plt.plot(r,MW.MassEnclosedTotal(r),label='Total Mass',linestyle='solid')
plt.xlabel('Radius (kpc)',fontsize='14')
plt.ylabel(r'Mass Enclosed ($M_{\odot}$)',fontsize='14')
plt.legend()
plt.title('Milky Way Mass Profile',fontsize='18')
plt.savefig('C:/Users/Jimmy/400B_Lilly/Homeworks/Homework5/MW Mass Profile.png')

# Plot simulated and Hernquist halo mass profiles of Milky Way
plt.figure(2)
plt.semilogy()
plt.plot(r,MW.MassEnclosed(1,r),label='Halo Mass ({0}'.format(str(MW_Mhalo))+r'$\cdot 10^{12} M_{\odot})$',linestyle='solid')
plt.plot(r,MW.HernquistM(r,MW_a,MW_Mhalo),label='Hernquist Mass (a={0})'.format(str(MW_a)),linestyle='dashed')
plt.xlabel('Radius (kpc)',fontsize='14')
plt.ylabel(r'Mass Enclosed ($M_{\odot}$)',fontsize='14')
plt.legend()
plt.title('Milky Way Hernquist - Halo Comparison',fontsize='18')
plt.savefig('C:/Users/Jimmy/400B_Lilly/Homeworks/Homework5/MW Hernq - Halo Comparison.png')   

# Plot rotation curve of Milky Way (each component and total)
plt.figure(3)
plt.plot(r,MW.CircularVelocity(1,r),label='Halo Velocity',linestyle='dashed')
plt.plot(r,MW.CircularVelocity(2,r),label='Disk Velocity',linestyle=':')
plt.plot(r,MW.CircularVelocity(3,r),label='Bulge Velocity',linestyle='dashdot')
plt.plot(r,MW.CircularVelocityTotal(r),label='Total Velocity',linestyle='solid')
plt.plot(r,MW.HernquistVCirc(r,MW_a,MW_Mhalo),label='Hernquist Velocity (a={0})'.format(str(MW_a)),linestyle='dashed')
plt.xlabel('Radius (kpc)',fontsize='14')
plt.ylabel('Circular Speed (km/s)',fontsize='14')
plt.legend()
plt.title('Milky Way Rotation Curve',fontsize='18')
plt.savefig('C:/Users/Jimmy/400B_Lilly/Homeworks/Homework5/MW Rotation Curve.png')
plt.show()
"""
"""
# Uncomment below to make plots for M31
# Call class for snapshot 0 of M31
M31 = MassProfile('M31',0)

# Initialize M31 scale radius (in kpc) and halo mass (from HW 3 in 10^12 Msun)
M31_a = 61.5*u.kpc
M31_Mhalo = 1.921

# Plot mass profile of M31 (each component and total)
plt.figure(1)
plt.semilogy()
plt.plot(r,M31.MassEnclosed(1,r),label='Halo Mass',linestyle='dashed')
plt.plot(r,M31.MassEnclosed(2,r),label='Disk Mass',linestyle=':')
plt.plot(r,M31.MassEnclosed(3,r),label='Bulge Mass',linestyle='dashdot')
plt.plot(r,M31.MassEnclosedTotal(r),label='Total Mass',linestyle='solid')
plt.xlabel('Radius (kpc)',fontsize='14')
plt.ylabel(r'Mass Enclosed ($M_{\odot}$)',fontsize='14')
plt.legend()
plt.title('M31 Mass Profile',fontsize='18')
plt.savefig('C:/Users/Jimmy/400B_Lilly/Homeworks/Homework5/M31 Mass Profile.png')

# Plot simulated and Hernquist halo mass profiles of M31
plt.figure(2)
plt.semilogy()
plt.plot(r,M31.MassEnclosed(1,r),label='Halo Mass ({0}'.format(str(M31_Mhalo))+r'$\cdot 10^{12} M_{\odot})$',linestyle='solid')
plt.plot(r,M31.HernquistM(r,M31_a,M31_Mhalo),label='Hernquist Mass (a={0})'.format(str(M31_a)),linestyle='dashed')
plt.xlabel('Radius (kpc)',fontsize='14')
plt.ylabel(r'Mass Enclosed ($M_{\odot}$)',fontsize='14')
plt.legend()
plt.title('M31 Hernquist - Halo Comparison',fontsize='18')
plt.savefig('C:/Users/Jimmy/400B_Lilly/Homeworks/Homework5/M31 Hernq - Halo Comparison.png')   

# Plot rotation curve of Milky Way (each component and total)
plt.figure(3)
plt.plot(r,M31.CircularVelocity(1,r),label='Halo Velocity',linestyle='dashed')
plt.plot(r,M31.CircularVelocity(2,r),label='Disk Velocity',linestyle=':')
plt.plot(r,M31.CircularVelocity(3,r),label='Bulge Velocity',linestyle='dashdot')
plt.plot(r,M31.CircularVelocityTotal(r),label='Total Velocity',linestyle='solid')
plt.plot(r,M31.HernquistVCirc(r,M31_a,M31_Mhalo),label='Hernquist Velocity (a={0})'.format(str(M31_a)),linestyle='dashed')
plt.xlabel('Radius (kpc)',fontsize='14')
plt.ylabel('Circular Speed (km/s)',fontsize='14')
plt.legend()
plt.title('M31 Rotation Curve',fontsize='18')
plt.savefig('C:/Users/Jimmy/400B_Lilly/Homeworks/Homework5/M31 Rotation Curve.png')
plt.show()
"""
"""
# Uncomment below to make plots for M33
# Call class for snapshot 0 of M33
M33 = MassProfile('M33',0)

# Initialize M33 scale radius (in kpc) and halo mass (from HW 3 in 10^12 Msun)
M33_a = 25*u.kpc
M33_Mhalo = .187

# Plot mass profile of M33 (each component and total)
plt.figure(1)
plt.semilogy()
plt.plot(r,M33.MassEnclosed(1,r),label='Halo Mass',linestyle='dashed')
plt.plot(r,M33.MassEnclosed(2,r),label='Disk Mass',linestyle=':')
plt.plot(r,M33.MassEnclosedTotal(r),label='Total Mass',linestyle='solid')
plt.xlabel('Radius (kpc)',fontsize='14')
plt.ylabel(r'Mass Enclosed ($M_{\odot}$)',fontsize='14')
plt.legend()
plt.title('M33 Mass Profile',fontsize='18')
plt.savefig('C:/Users/Jimmy/400B_Lilly/Homeworks/Homework5/M33 Mass Profile.png')

# Plot simulated and Hernquist halo mass profiles of M33
plt.figure(2)
plt.semilogy()
plt.plot(r,M33.MassEnclosed(1,r),label='Halo Mass ({0}'.format(str(M33_Mhalo))+r'$\cdot 10^{12} M_{\odot})$',linestyle='solid')
plt.plot(r,M33.HernquistM(r,M33_a,M33_Mhalo),label='Hernquist Mass (a={0})'.format(str(M33_a)),linestyle='dashed')
plt.xlabel('Radius (kpc)',fontsize='14')
plt.ylabel(r'Mass Enclosed ($M_{\odot}$)',fontsize='14')
plt.legend()
plt.title('M33 Hernquist - Halo Comparison',fontsize='18')
plt.savefig('C:/Users/Jimmy/400B_Lilly/Homeworks/Homework5/M33 Hernq - Halo Comparison.png')   

# Plot rotation curve of M33 (each component and total)
plt.figure(3)
plt.plot(r,M33.CircularVelocity(1,r),label='Halo Velocity',linestyle='dashed')
plt.plot(r,M33.CircularVelocity(2,r),label='Disk Velocity',linestyle=':')
plt.plot(r,M33.CircularVelocityTotal(r),label='Total Velocity',linestyle='solid')
plt.plot(r,M33.HernquistVCirc(r,M33_a,M33_Mhalo),label='Hernquist Velocity (a={0})'.format(str(M33_a)),linestyle='dashed')
plt.xlabel('Radius (kpc)',fontsize='14')
plt.ylabel('Circular Speed (km/s)',fontsize='14')
plt.legend()
plt.title('M33 Rotation Curve',fontsize='18')
plt.savefig('C:/Users/Jimmy/400B_Lilly/Homeworks/Homework5/M33 Rotation Curve.png')
plt.show()
"""