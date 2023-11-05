import numpy as np
from sympy import *

n_consumers = 5
n_resources = 2
alpha_mat = np.array([[Symbol("alpha_{}_{}".format(i,j)) for j in range(n_resources)] for i in range(n_resources)])
beta_mat = np.array([[Symbol("beta_{}_{}".format(i,j)) for j in range(n_resources)] for i in range(n_consumers)])
C_vec = np.array([Symbol("C_{}".format(i)) for i in range(n_consumers)])
R_vec = np.array([Symbol("R_{}".format(i)) for i in range(n_resources)])
m_vec = np.array([Symbol("m_{}".format(i)) for i in range(n_consumers)])
p_vec = np.array([Symbol("p_{}".format(i)) for i in range(n_consumers)])
b = Symbol('b')

omega = 0
for i in range(n_consumers):
    for j in range(n_resources):
        omega += beta_mat[i][j]*R_vec[j]

print(omega)

f_vec = np.array([p_vec[i]*C_vec[i]] for i in range(n_consumers))
print(f_vec)

print(alpha_mat)
print(beta_mat)
print(C_vec)
print(R_vec)