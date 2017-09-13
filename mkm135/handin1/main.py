import numpy as np
import matplotlib.pyplot as plt


N_epi = 1e15
L_n_epi = 400e-6
N_sub = np.array([1e15, 1e16, 1e17], dtype=np.float64)
t_p = 200e-9
t_epi = 4e-6
t_sub = 250e-6
T = 300     # temperature, K
A = 0.1e-6  # cross-section area, m^2
k = 1.38e-23
q = 1.6e-19
V_T = k*T/q

L_n_sub = L_n_epi*N_sub/N_epi
N_p = (N_epi*t_epi + N_sub*t_sub) / t_p
L_p = L_n_epi*N_p/N_epi

N_A = N_p
N_D = N_sub + N_epi
L_P = L_p
L_N = L_n_sub + L_n_epi

n_i = 10e10  # lecture 1&2 slide 20
mu_n = np.array([1.5e3, 1.3e3, 8e2], dtype=np.float64)
mu_p = np.array([3.8e2, 3.5e2, 2e2], dtype=np.float64)

D_N = V_T * mu_n
D_P = V_T * mu_p

V = 2  # np.linspace(0, 2, 101)

for i in range(3):
    plt.plot(V, q*(D_N/L_N * n_i**2/N_A + D_P/L_P * n_i**2/N_D)*(np.exp(V*V_T) - 1))
plt.show()
