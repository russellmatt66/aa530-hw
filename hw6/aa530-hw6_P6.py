"""
Matt Russell
University of Washington
Department of Aeronautics & Astronautics
11/30/21
AA530: Mechanics of Solids
HW6 Problem 6
Plot Von Mises Yield surface
Explicit expression for sigma_{3} = sigma_{3}(sigma_{1},sigma_{2},\tau_{i}) is obtainable
Second-order polynomial, solvable with quadratic formula
"""
import numpy as np
import matplotlib.pyplot as plt

def f_3hat(sig1_hat,sig2_hat,porm): # expression for normalized principal stress (in 'z'-dir)
    # Solution is based on quadratic formula, therefore 'porm' variable indicates whether to take positive or negative solution
    if porm == 'p': # positive solution
        return sig1_hat + sig2_hat + np.sqrt(-np.power(sig1_hat,2)-np.power(sig2_hat,2) + 4.0*np.multiply(sig1_hat,sig2_hat) + 2)
    if porm == 'n': # negative solution
        return sig1_hat + sig2_hat - np.sqrt(-np.power(sig1_hat,2)-np.power(sig2_hat,2) + 4.0*np.multiply(sig1_hat,sig2_hat) + 2)
    else:
        print('"porm" argument should be either "p" or "n"')
        return None

sigma_y = 200.0 # [MPa] Yield stress

nsigma = 100

sig1_hat = np.linspace(-1.0,1.0,nsigma) 
sig2_hat = np.linspace(-1.0,1.0,nsigma)

ss1_hat, ss2_hat = np.meshgrid(sig1_hat,sig2_hat)

ss3_hat_p = f_3hat(ss1_hat,ss2_hat,'p')
ss3_hat_n = f_3hat(ss1_hat,ss2_hat,'n')

negtest_fig = plt.figure() 
negtest_ax = plt.axes(projection='3d')
negtest_ax.plot_surface(ss1_hat,ss2_hat,ss3_hat_p)
negtest_ax.plot_surface(ss1_hat,ss2_hat,ss3_hat_n)
negtest_ax.set_title('Von Mises Yield')
negtest_ax.set_xlabel('$\\tilde{\sigma}_{1}$')
negtest_ax.set_ylabel('$\\tilde{\sigma}_{2}$')
negtest_ax.set_zlabel('$\\tilde{\sigma}_{3}$')
negtest_ax.set_xlim3d([-1.0,1.0])
negtest_ax.set_ylim3d([-1.0,1.0])

"""
#negtest_ax.contour3D(ss1_hat,ss2_hat,ss3_hat_n)
negtest_ax.plot_surface(ss1_hat,ss2_hat,ss3_hat_n)
negtest_ax.set_title('Negative solution')
"""
"""
negtest_ax.contour3D(ss1_hat,ss2_hat,ss3_hat_p)
negtest_ax.set_title('Positive solution')
"""
plt.show()
