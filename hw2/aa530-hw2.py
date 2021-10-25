"""
Matt Russell
AA530: Mechanics of Solids
HW2
10/24/21
"""
import numpy as np
import matplotlib.pyplot as plt

""" Problem 1 """
# 1.3 and 1.4
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

# 1.6
alpha = np.linspace(0.0,180.0,100) 
alpharad = np.deg2rad(alpha) # angle in radians for convenience
Tn = (1.0/2.0)*(sigma_1 + sigma_2) + 0.5*(sigma_1 - sigma_2)*np.cos(2.0*alpharad) + 0.5*(tau_12 + tau_21)*np.sin(2.0*alpharad)
Ttau = 0.5*(sigma_1 - sigma_2)*np.sin(2.0*alpharad) - 0.5*(tau_12 + tau_21)*np.cos(2.0*alpharad)

fig_stress1pt6 = plt.figure()
plt.plot(Tn,Ttau)
plt.title('Stress Trajectory')
plt.xlabel('Normal Stress [MPa]')
plt.ylabel('Shear Stress [MPa]')

plt.show()
