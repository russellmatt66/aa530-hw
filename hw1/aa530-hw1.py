"""
Matt Russell
AA 530: Solid Mechanics
10/11/21
"""
import numpy as np

""" Problem 2.5 """
LST_2pt5 = np.array([[0.0,0.0,0.0],[0.0,0.0,0.25],[0.0,0.25,0.125]])
evals, evecs = np.linalg.eig(LST_2pt5)
print("The eigenvalues of the LST in 2.5 are %f", evals)
print("The eigenvectors of the LST in 2.5 are %f", evecs)