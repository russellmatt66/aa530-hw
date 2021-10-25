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
sigma_avg = 0.5*(sigma_1 + sigma_2) # [MPa]
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
plt.plot(Tn,Ttau,label='Stress')
plt.scatter(sigma_avg,0,label='Center')
plt.xlim([0,35])
plt.ylim([-16,16])
plt.title('Stress Trajectory')
plt.xlabel('Normal Stress [MPa]')
plt.ylabel('Shear Stress [MPa]')
plt.legend()

# 1.7
T_nprin1 = sigma_avg + np.sqrt(0.25*(sigma_1 - sigma_2)**2 + tau_12**2)
T_nprin2 = sigma_avg - np.sqrt(0.25*(sigma_1 - sigma_2)**2 + tau_12**2)
T_taumaxpos = np.sqrt(0.25*(sigma_1 - sigma_2)**2 + tau_12**2)
T_taumaxneg = -np.sqrt(0.25*(sigma_1 - sigma_2)**2 + tau_12**2)

print("The principal stresses for 1.6 are %f and %f [MPa]" %(T_nprin1,T_nprin2))
print("The maximum shear stress for Problem 1 is %f [MPa]" %T_taumaxpos) 

# 2.2 - principal stresses and directions of stress tensor
sigmatensor_2pt2 = np.array([[6.0, -2.0, 0.0],[-2.0, 3.0, 4.0],[0.0,4.0,3.0]])
evals_2pt2,evecs_2pt2 = np.linalg.eig(sigmatensor_2pt2)

print("The eigenvalues (principal stresses) are", evals_2pt2)
print("The eigenvectors (principal directions) are\n", evecs_2pt2)

# 2.3 - maximum shear stress
T_taumax2pt3 = np.sqrt(0.25*(evals_2pt2[0] - evals_2pt2[1])**2 + 4) # 4 - tau_{xy}^{2}
print("The maximum shear stress is %f [Pa]" %T_taumax2pt3)

# 2.4 - Von Mises stress
sigma_VM = np.sqrt(0.5*((evals_2pt2[0] - evals_2pt2[1])**2 + (evals_2pt2[0] - evals_2pt2[2])**2 + (evals_2pt2[1] - evals_2pt2[2])**2))
print("The Von Mises stress is %f [Pa]" %sigma_VM)
    
plt.show()
