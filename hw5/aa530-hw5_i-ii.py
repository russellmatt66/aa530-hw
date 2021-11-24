"""
Matt Russell
University of Washington
Department of Aeronautics & Astronautics
AA530 HW5
Problems 1 and 2
Calculate cauchy and nominal stress in principal (e1) direction as a function of stretch rate for uniaxial and biaxial strain states
"""
import numpy as np
import matplotlib.pyplot as plt

# Representative Material Properties: Table 3.7 in Allan and Bower
mu_1_nh = 0.4 # MN m^{-2}, Neo-Hookean

mu_1_mr = 0.39 # MN m^{-2}, Mooney-Rivlin
mu_2_mr = 0.015 # MN m^{-2}, Mooney-Rivlin

mu_1_ab = 0.4 # MN m^{-2}, Arruda-Boyce
beta_ab = 10.0 # Arruda-Boyce material parameter

mu_og = np.array([0.62,0.00118,-0.00981]) # MN m^{-2}, Ogden
alpha_og = np.array([1.3,5.0,-2.0]) # Ogden material parameters

# Stretch Rate
nl = 100
lambda_sr = np.linspace(0.5,2.0,nl)

""" Problem 1 - Uniaxial tension """
I_1_uni = lambda_sr**2 + 2.0/lambda_sr

# Cauchy

# Nominal - Table 3.6 in Allan and Bower
S1_nh_uni = mu_1_nh*(lambda_sr - 1.0/lambda_sr**2)

S1_mr_uni = mu_1_mr*(lambda_sr - 1.0/lambda_sr**2) + mu_2_mr*(1.0 - 1.0/lambda_sr**3)

C_ab_uni = mu_1_ab*(1.0 + I_1_uni/(5.0*beta_ab**2) + (33.0*I_1_uni**2)/(525.0*beta_ab**4))
S1_ab_uni = C_ab_uni*(lambda_sr - 1.0/lambda_sr**2)

S1_og_uni = 0.0
for oidx in np.arange(mu_og.shape[0]):
    S1_og_uni += mu_og[oidx]*(lambda_sr**(alpha_og[oidx]) - 1.0/lambda_sr**(-alpha_og[oidx]/2)) * (1.0/lambda_sr)

""" Problem 2 - Biaxial tension """
I_1_bi = 2.0*lambda_sr**2 + 1.0/lambda_sr**4

# Cauchy

# Nominal - Table 3.6 in Allan and Bower
S1_nh_bi = mu_1_nh*(lambda_sr - 1.0/lambda_sr**5)

S1_mr_bi = mu_1_mr*(lambda_sr - 1.0/lambda_sr**5) + mu_2_mr*(lambda_sr**3 - 1.0/lambda_sr**3)

C_ab_bi = mu_1_ab*(1.0 + I_1_bi/(5.0*beta_ab**2) + (33.0*I_1_bi**2)/(525.0*beta_ab**4))
S1_ab_bi = C_ab_bi*(lambda_sr - 1.0/lambda_sr**5)

S1_og_bi = 0.0
for oidx in np.arange(mu_og.shape[0]):
    S1_og_bi += mu_og[oidx]*(lambda_sr**(alpha_og[oidx]) - lambda_sr**(-2.0*alpha_og[oidx])) * (1.0/lambda_sr)

""" Graphs """
# Problem 1

# Problem 2
