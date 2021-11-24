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

# Deformation characteristics + tools (uniform deformation), i.e, tensor stuff
F = np.eye(3,3) # Deformation-gradient 
B = np.eye(3,3) # Left Cauchy-Green
J = 1.0 # Jacobian, J = det(F), implies that invariants are already normalized

BddB = 0.0 # shows up in I_2, pg. 96 of Allan and Bower
for idx in np.arange(B.shape[0]):
    for kdx in np.arange(B.shape[1]):
        BddB = BddB + B[idx,kdx]*B[kdx,idx]


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
# Invariants, J = 1 so already normalized
I_1_uni = lambda_sr**2 + 2.0/lambda_sr 
I_2_uni = 0.5*(I_1_uni**2 - BddB)

# Cauchy
sigma11_nh_uni = mu_1_nh/J**(5/3)*(B[0,0] - (1.0/3.0)*np.trace(B)) # Neo-Hookean

## Mooney-Rivlin
# tensor contraction
B1k_Bk1 = 0.0
for kidx in np.arange(B.shape[0]): # B is square
    B1k_Bk1 = B1k_Bk1 + B[0,kidx]*B[kidx,0]
    
sigma11_mr_uni = (mu_1_mr/J**(5/3)) * (1.0 + (1.0/J**(2/3))*np.trace(B)) * B[0,0] \
    - (mu_2_mr/(3.0*J**(5/3))) * (I_1_uni + 2.0*I_2_uni) - (mu_2_mr/J**(7/3)) * B1k_Bk1

## Arruda-Boyce
sigma11_ab_uni = (2.0/J**(5/3))*(0.5*mu_1_ab + (1.0/(10.0*beta_ab**2))*I_1_uni) * B[0,0] \
    - (2.0/(3.0*J))*(I_1_uni * (0.5*mu_1_ab + (1.0/(10.0 * beta_ab**2))*I_1_uni))

# Ogden defeated me for the moment
"""
pdvSED = 0.0 # Eq (31) in HW4
for kidx in np.arange(mu_og.shape[0]):
    pdvSED = pdvSED + (2.0*mu_og[kidx])/(alpha_og[kidx]**2)*(1.0/J**(1/3))

# sigma11_og_uni = (1.0/J**(1/3))# Equation (32) in HW4, outer product (which has a typo) should just contribute 1 to principal component
"""

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
I_2_bi = 0.5*(I_1_bi**2 - BddB)

# Cauchy
sigma11_nh_bi = mu_1_nh/J**(5/3)*(B[0,0] - (1.0/3.0)*np.trace(B)) # Neo-Hookean

## Mooney-Rivlin
sigma11_mr_bi = (mu_1_mr/J**(5/3)) * (1.0 + (1.0/J**(2/3))*np.trace(B)) * B[0,0] \
    - (mu_2_mr/(3.0*J**(5/3))) * (I_1_bi + 2.0*I_2_bi) - (mu_2_mr/J**(7/3)) * B1k_Bk1

## Arruda-Boyce
sigma11_ab_bi = (2.0/J**(5/3))*(0.5*mu_1_ab + (1.0/(10.0*beta_ab**2))*I_1_bi) * B[0,0] - \
    (2.0/(3.0*J))*(I_1_bi * (0.5*mu_1_ab + (1.0/(10.0 * beta_ab**2))*I_1_bi))

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
## Subplots
nr = 2
nc = 2
nh_fig,nh_ax = plt.subplots(nr,nc)
mr_fig,mr_ax = plt.subplots(nr,nc)
ab_fig,ab_ax = plt.subplots(nr,nc)

sigma11_nh_uni = sigma11_nh_uni * np.ones(lambda_sr.shape[0])
sigma11_nh_bi = sigma11_nh_bi * np.ones(lambda_sr.shape[0])

# Neo - Hookean
nh_ax[0,0].plot(lambda_sr,sigma11_nh_uni,label='i')
nh_ax[0,1].plot(lambda_sr,S1_nh_uni,label='ii')
nh_ax[1,0].plot(lambda_sr,sigma11_nh_bi,label='iii')
nh_ax[1,1].plot(lambda_sr,S1_nh_bi,label='iv')

nh_ax[0,0].legend()
nh_ax[0,1].legend()
nh_ax[1,0].legend()
nh_ax[1,1].legend()

# Mooney-Rivlin
nh_ax[0,0].plot(lambda_sr,sigma11_nh_uni,label='i')
nh_ax[0,1].plot(lambda_sr,S1_nh_uni,label='ii')
nh_ax[1,0].plot(lambda_sr,sigma11_nh_bi,label='iii')
nh_ax[1,1].plot(lambda_sr,S1_nh_bi,label='iv')

nh_ax[0,0].legend()
nh_ax[0,1].legend()
nh_ax[1,0].legend()
nh_ax[1,1].legend()

# Arruda-Boyce
nh_ax[0,0].plot(lambda_sr,sigma11_nh_uni,label='i')
nh_ax[0,1].plot(lambda_sr,S1_nh_uni,label='ii')
nh_ax[1,0].plot(lambda_sr,sigma11_nh_bi,label='iii')
nh_ax[1,1].plot(lambda_sr,S1_nh_bi,label='iv')

nh_ax[0,0].legend()
nh_ax[0,1].legend()
nh_ax[1,0].legend()
nh_ax[1,1].legend()
plt.show()
