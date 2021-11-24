"""
Matt Russell
University of Washington
Department of Aeronautics & Astronautics
AA530 HW5 
Problem 4
"""
import numpy as np
import matplotlib.pyplot as plt

sigma_0 = 1.0 # stress load [MPa]
k = 0.1e-3 # elastic constant [MPa]
eta = 20.0e-3 # plastic constant [MPa s]

nt = 100
t = np.linspace(0.0,600.0,nt)

epsilon = sigma_0*(1.0/k + 1.0/eta * t) # strain []

plt.plot(t,epsilon)
plt.xlabel('time (s)')
plt.ylabel('epsilon')
plt.xlim((t[0],t[nt-1]))
plt.ylim((0.0,epsilon[nt-1]))
plt.title('Strain of a viscoelastic material with k = 0.1 GPa, $\eta$ = 20 GPa$\cdot$s, $\sigma_{0}$ = 1 MPa') 

plt.show()
