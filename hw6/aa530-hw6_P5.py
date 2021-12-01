"""
Matt Russell
University of Washington
Department of Aeronautics & Astronautics
AA530: Mechanics of Solids HW6
11/29/21
Plotting of Tresca yield surface 
O(N^{3}) solution for the yield criteria - doesn't work, 3D coordinate array is wrong size for meshgrid()
Using analytical expression from notes instead
"""
import numpy as np
import matplotlib.pyplot as plt

# Problem 5 - Plot Tresca yield surface
sigma_yield = 200.0 # [MPa] uniaxial yield stress


"""
O(N^3) solution for posterity

def maximum(a,b):
    if a >= b:
        return a
    else:
        return b
        
nsigma = 25

sigma_pr1 = np.linspace(0.0,sigma_yield,nsigma) # vectorized principal stresses
sigma_pr2 = np.linspace(0.0,sigma_yield,nsigma)
sigma_pr3 = np.linspace(0.0,sigma_yield,nsigma)

sigma_surf1 = np.empty(1) # Containers for yield surface points 
sigma_surf2 = np.empty(1)
sigma_surf3 = np.empty(1)

tol = 1.0e-13
for idx1 in np.arange(sigma_pr1.shape[0]): # O(N^3)
    print('Loop %i' %idx1)
    for idx2 in np.arange(sigma_pr2.shape[0]):
        for idx3 in np.arange(sigma_pr3.shape[0]):
            max_intrmd = maximum(np.abs(sigma_pr1[idx1] - sigma_pr2[idx2]),np.abs(sigma_pr1[idx1] - sigma_pr3[idx3]))
            max_real = maximum(max_intrmd,np.abs(sigma_pr2[idx2] - sigma_pr3[idx3]))
            #max_intrmd = np.maximum(np.abs(sigma_pr1[idx1] - sigma_pr2[idx2]),np.abs(sigma_pr1[idx1] - sigma_pr3[idx3]))
            #max_real = np.maximum(max_intrmd,np.abs(sigma_pr2[idx2] - sigma_pr3[idx3]))
            if np.abs(max_real - sigma_yield) < tol:
                sigma_surf1 = np.append(sigma_surf1,sigma_pr1[idx1])
                sigma_surf2 = np.append(sigma_surf2,sigma_pr2[idx2])
                sigma_surf3 = np.append(sigma_surf3,sigma_pr3[idx3])
 
sigma_surf1 = np.delete(sigma_surf1,0)
sigma_surf2 = np.delete(sigma_surf2,0)
sigma_surf3 = np.delete(sigma_surf3,0)

print(sigma_surf1.shape)
print(sigma_surf2.shape)
print(sigma_surf3.shape)

spr1, spr2, spr3 = np.meshgrid(sigma_surf1, sigma_surf2, sigma_surf3) # coordinate arrays for the surface points

print(spr1.shape)
print(spr2.shape)
print(spr3.shape)

P5_fig = plt.figure()
P5_ax = plt.axes(projection='3d')
P5_ax.plot_surface(spr1,spr2,spr3) # Doesn't work b/c arrays are 3D, use analytical expression instead

plt.show()
 