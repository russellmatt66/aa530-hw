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
d = 3
gamma = 0.5
F = np.eye(d) # Deformation Gradient Tensor
F[0,1] = gamma

print("The Deformation Gradient Tensor, F, is\n")
print(F)

J = np.linalg.det(F) # Jacobian

print("The Jacobian of F is %4.3f" %J)

B = np.matmul(F,np.transpose(F)) # Left Cauchy-Green
"""
There's a bug in the calculation of B from principal stretches and directions but both ways give the same B so it must be downstream 
nr = 3
nc = 3
B = np.zeros((nr,nc))
for i in np.arange(nr):
    for j in np.arange(nc):
        B[i,j] = np.dot(F[i,:],F[j,:])
"""
C = np.matmul(np.transpose(F),F) # Right "

print("1.1. The Left Cauchy-Green tensor of F, B, is\n")
print(B)
print("1.1. The Right Cauchy-Green tensor of F, C, is\n")
print(C)

""" Problem 1.2 - Eigenvalues of B and C """
evals_B, evecs_B = np.linalg.eig(B)
evals_C, evecs_C = np.linalg.eig(C)

print("1.2. The eigenvalues of B are\n")
print(evals_B)
print("1.2. The eigenvalues of C are\n")
print(evals_C)

""" Problem 1.3 - Principal Stretches and Directions of B """
lambda_B = np.sqrt(evals_B) # {\lambda_{i} | \lambda_{i} = sqrt(e_{i})}

print("1.3. The principal stretches (scaled eigenvalues) of B are\n")
print(lambda_B)
print("1.3. The principal directions, i.e, eigenvectors of B are\n")
print(evecs_B)

""" Problem 1.4 - Verification """
B_verif = np.zeros((B.shape[0],B.shape[1]))

for eidx in np.arange(np.size(evals_B)):
    B_verif = B_verif + lambda_B[eidx]**2 * np.outer(evecs_B[eidx],evecs_B[eidx])

print("1.4. The Left Cauchy-Green tensor of F, B, is\n")
print(B)
print("1.4. The Left Cauchy-Green tensor of F, B, that was constructed from the principal stretches and directions is\n")
print(B_verif)

""" Problem 1.5 - Invariants and Normalized Forms """
# Invariants
I_1 = np.trace(B)
I_2 = 0.5*(I_1**2 - np.trace(np.matmul(B,B)))
I_3 = np.linalg.det(B)
I_array = np.asarray([I_1, I_2, I_3])
for iidx in np.arange(I_array.shape[0]):
    print("Invariant %i is %f" %(iidx,I_array[iidx]))

# Normalized Forms
Ibar_1 = I_1 / J**(2/3)
Ibar_2 = I_2 / J**(4/3)
Ibar_3 = 1
Ibar_array = np.asarray([Ibar_1, Ibar_2, Ibar_3])
for bidx in np.arange(Ibar_array.shape[0]):
    print("Normalized invariant %i is %f" %(bidx,Ibar_array[bidx]))

