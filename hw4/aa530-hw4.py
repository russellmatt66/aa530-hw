"""
Matt Russell
University of Washington
Department of Aeronautics & Astronautics
AA530: Mechanics of Solids
HW4
11/8/21
"""
import numpy as np

""" Problem 1.1 - Left and Right Cauchy-Green Tensors """
gamma = 0.5
F = np.eye((3,3))
F[0,1] = gamma

B = np.matmul(F,F.transpose)
C = np.matmul(F.transpose,F)

print("1.1. The Left Cauchy-Green tensor of F, B, is %f" %B)
print("1.1. The Right Cauchy-Green tensor of F, C, is %f" %C)

""" Problem 1.2 - Eigenvalues of B and C """
evals_B, evecs_B = np.linalg.eig(B)
evals_C, evecs_C = np.linalg.eig(C)

print("1.2. The eigenvalues of B are %f" %evals_B)
print("1.2. The eigenvalues of C are %f" %evals_C)

""" Problem 1.3 - Principal Stretches and Directions of B """
prncplstrtch_B = np.sqrt(evals_B) # {\lambda_{i} | \lambda_{i} = sqrt(e_{i})}

print("1.3. The principal stretches of B are %f" %prncplstrtch_B)
print("1.3. The principal directions of B are %f" %evecs_B)

""" Problem 1.4 - Verification """
B_verif = np.zeros((B.shape[0],B.shape[1]))
for eidx in np.size(evals_B):
    B_verif = B_verif + evals_B[eidx] * np.outer(evecs_B[eidx],evecs_B[eidx])

print("1.4. The Left Cauchy-Green tensor of F, B, is %f" %B)
print("1.4. The Left Cauchy-Green tensor of F, B, constructed from the principal stretches and directions is %f" %B_verif)

""" Problem 1.5 - Invariants and Normalized Form """
