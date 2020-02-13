# Program to test my COM component calculation

# Import nump
import numpy as np

# Define x positions of particles and their masses 
# (one particle has large mass to verify COM component should be closest to it)
x = [1,2,3]
mass = [4,5,10000]

# Find product of x-component of positions and their respective masses
xmProduct = np.multiply(x,mass)

# Calculate the x-component of the center of mass
XCOM = np.sum(xmProduct)/np.sum(mass)

print(XCOM)