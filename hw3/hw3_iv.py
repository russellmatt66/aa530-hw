"""
Matt Russell
AA530 HW3 - Problem 4
Dept of Aero&Astro, UW
11/1/21

Compute stresses based on readings from a strain rosette
"""
import numpy as np

E = 200.0 # [GPa]
nu = 0.3

epsilon_x = 0.005
epsilon_y = -0.001
epsilon_xy = 3.0/1732.0

sigma_x = (E)/((1.0 + nu)*(1.0 - 2.0*nu))*((1.0-nu)*epsilon_x + nu*epsilon_y)
sigma_y = (E)/((1.0 + nu)*(1.0 - 2.0*nu))*(nu*epsilon_x + (1.0 - nu)*epsilon_y)
sigma_xy = (E)/((1.0 + nu)*(1.0 - 2.0*nu))*((1.0 - 2.0*nu)/(2.0))*2.0*epsilon_xy

print("sigma_x is %f" %sigma_x)
print("sigma_y is %f" %sigma_y)
print("sigma_xy is %f" %sigma_xy)
