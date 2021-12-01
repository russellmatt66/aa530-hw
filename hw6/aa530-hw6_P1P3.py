"""
Matt Russell
University of Washington
Department of Aeronautics & Astronautics
AA530: Mechanics of Solids HW6
Simple script plotting direct relationships
11/29/21
"""
import numpy as np
import matplotlib.pyplot as plt

# Problem 1 - Stress in a viscoelastic material (Maxwell)
eps_0 = 0.1 # Strain
k = 0.1 # [GPa]
eta = 20.0 # [GPa-s]

nt1 = 1000
t1 = np.linspace(0.0,600.0,nt1)

sigma_1 = k * eps_0 * np.exp(-(k/eta)*t1) # Constitutive relationship from HW5.5

P1_fig = plt.figure()
plt.plot(t1,sigma_1)
plt.title('Viscoelastic material (Maxwell model) subject to a constant strain')
plt.xlabel('time (s)')
plt.ylabel('Stress (GPa)')

# Problem 3 - Strain in a viscoelastic material (Kelvin-Voigt)
sigma_0 = 1.0e-3 # [MPa] -> [GPa]

eps_3 = (sigma_0/k) * (1.0 - np.exp(-(k/eta) * t1))

P3_fig = plt.figure()
plt.plot(t1,eps_3)
plt.title('Viscoelastic material (Kelvin-Voigt model) subject to a constant stress')
plt.xlabel('time (s)')
plt.ylabel('Strain')

plt.show()
 