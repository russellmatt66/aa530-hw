"""
Matt Russell
AA530 HW3
Department of Aero&Astro, UW
11/1/21
"""
import numpy as np

def isothermalConstvReln(epsilon,E,nu):
    # epsilon - strain state at a point in a material, 3x3 matrix (rank-2 tensor)
    # Function returns the corresponding stress tensor at a point in an isotropic material
    # Performs transformation in O(N^2) (N=3)
    sigma = np.empty((np.size(epsilon,axis=0),np.size(epsilon,axis=1)))
    for i in np.arange(sigma.shape[0]):
        for j in np.arange(sigma.shape[1]):
            sigma[i,j] = (E/(1.0 + nu))*epsilon[i,j]
            if i == j:
                sigma[i,j] += (nu/(1.0 - 2.0*nu))*np.trace(epsilon)
    return sigma

epsilon_matrixform = np.asarray([[0.001,0.001,-0.002],[0.001,0.0,0.005],[-0.002,0.005,0.0]]) # state of strain at a point in an isotropic material

youngModulus = 200.0 # [GPa]
poiRat = 0.3 # Poisson's ratio

sigma = isothermalConstvReln(epsilon_matrixform,youngModulus,poiRat)
print(sigma)
                                         
