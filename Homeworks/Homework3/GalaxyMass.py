# Program to calculate the total mass of any galaxy component 
# (1=Halo,2=Disk,3=Bulge)

# Import important packages
import numpy as np
from ReadFile import Read
import astropy.units as u

def ComponentMass(file,partType):
    # Inputs
    #    file = text file to extract data from (galaxy and time of snapshot)
    #    partType = type of particle (see above)
    # Returns:
    #    final_component_mass = total mass of any galaxy component
    
    # Call 'Read' function from ReadFile.py
    time, numParticle, data = Read(file)
    
    # Get index for all particles of given type
    index = np.where(data['type']==partType)

    # Store particle masses of given type in an array (units of 10^12 Msun)
    component_mass = data['m'][index]/(10**2)*u.solMass
    
    # Calculate (and round) sum of all particle masses of given type
    final_component_mass = np.round(np.sum(component_mass),3)

    # Conditionals to recognize which galaxy we're investigating
    # Used in print statement below
    if 'MW' in file:
        string1 = 'Milky Way'
    elif 'M31' in file:
        string1 = 'M31'
    else:
        string1 = 'M33'

    # Conditionals to recognize which galaxy component mass we're finding
    # Used in print statement below
    if partType == 1:
        string2 = 'Halo'
    elif partType == 2:
        string2 = 'Disk'
    else:
        string2 = 'Bulge'

    # Print mass of galaxy component
    #print('Mass of {0} {1} = {2}'.format(string1,string2,final_component_mass))
    
    return final_component_mass

# The below code can be to used to calculate galaxy component masses, total
# masses, and baryon fractions.

# Use ComponentMass to compute the mass of the MW components
MW_halo_mass = ComponentMass('MW_000.txt',1)
MW_disk_mass = ComponentMass('MW_000.txt',2)
MW_bulge_mass = ComponentMass('MW_000.txt',3)

# Find MW stellar (bulk+disk) & total (halo+bulk+disk) masses & baryon fraction
MW_stellar_mass = np.round(MW_disk_mass+MW_bulge_mass,3)
MW_total_mass = np.round(MW_halo_mass+MW_stellar_mass,3)
MW_baryon_fraction = np.round(MW_stellar_mass/MW_total_mass,3)*u.solMass

# Use ComponentMass to compute the mass of the MW components
M31_halo_mass = ComponentMass('M31_000.txt',1)
M31_disk_mass = ComponentMass('M31_000.txt',2)
M31_bulge_mass = ComponentMass('M31_000.txt',3)

# Find M31 stellar (bulk+disk) & total (halo+bulk+disk) masses & baryon fraction
M31_stellar_mass = np.round(M31_disk_mass+M31_bulge_mass,3)
M31_total_mass = np.round(M31_halo_mass+M31_stellar_mass,3)
M31_baryon_fraction = np.round(M31_stellar_mass/M31_total_mass,3)*u.solMass

# Use ComponentMass to compute the mass of the MW components
M33_halo_mass = ComponentMass('M33_000.txt',1)
M33_disk_mass = ComponentMass('M33_000.txt',2)
M33_bulge_mass = ComponentMass('M33_000.txt',3)

# Find M33 stellar (bulk+disk) & total (halo+bulk+disk) masses & baryon fraction
M33_stellar_mass = np.round(M33_disk_mass+M33_bulge_mass,3)
M33_total_mass = np.round(M33_halo_mass+M33_stellar_mass,3)
M33_baryon_fraction = np.round(M33_stellar_mass/M33_total_mass,3)*u.solMass

# Calculate component & total masses for Local Group and its baryon fraction
LocalGroup_halo_mass = np.round(MW_halo_mass + M31_halo_mass + M33_halo_mass,3)
LocalGroup_disk_mass = np.round(MW_disk_mass + M31_disk_mass + M33_disk_mass,3)
LocalGroup_bulge_mass = np.round(MW_bulge_mass + M31_bulge_mass + M33_bulge_mass,3)
LocalGroup_mass = np.round(MW_total_mass+M31_total_mass+M33_total_mass,3)
LocalGroup_baryon_fraction = np.round((MW_stellar_mass+M31_stellar_mass+M33_stellar_mass)/LocalGroup_mass,3)

# Print final results (all in units of 10^12 Msun)
print(r'All results in units of 10^12 solar masses')
print('Total MW Mass = ' + str(MW_total_mass))
print('Total M31 Mass = ' + str(M31_total_mass))
print('Total M33 Mass = ' + str(M33_total_mass))
print('Total Local Group Mass = ' + str(LocalGroup_mass))

print('MW Baryon Fraction = ' + str(MW_baryon_fraction))
print('M31 Baryon Fraction = ' + str(M31_baryon_fraction))
print('M33 Baryon Fraction = ' + str(M33_baryon_fraction))
print('Local Group Baryon Fraction = ' + str(LocalGroup_baryon_fraction))

print('Local Group Halo Mass = ' + str(LocalGroup_halo_mass))
print('Local Group Disk Mass = ' + str(LocalGroup_disk_mass))
print('Local Group Bulge Mass = ' + str(LocalGroup_bulge_mass))