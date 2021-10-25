"""
Matt Russell
AA530: Mechanics of Solids
HW2
10/24/21
"""
import numpy as np
import matplotlib.pyplot as plt

""" Problem 1 """
# 1.3 
alpha = np.linspace(1.0,179.0,100) # degrees
sigma_1 = 30.0 # [MPa]
sigma_2 = 10.0 # [MPa]
tau_12 = -10.0 # [MPa]
tau_21 = tau_12 # [MPa]

fig_cot1pt3 = plt.figure()
plt.plot(alpha,np.cos(2.0*np.deg2rad(alpha))/np.sin(2.0*np.deg2rad(alpha)),label='cot(2$\\alpha$)')
plt.plot(alpha,np.tan(2.0*np.deg2rad(alpha)),label='tan(2$\\alpha$)')
plt.title('Problem 1.3 and 1.4')
plt.xlim([0.0,180.0])
plt.xlabel('$\\alpha$ (deg)$')
plt.legend()

plt.show()
