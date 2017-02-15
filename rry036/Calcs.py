import numpy as np
import matplotlib.pyplot as plt

f = 4e4
w = 2*np.pi*f
sigma = 4
mu_0 = 4*np.pi*1e-7
e_0 = 8.85e-12
e_r = 80
e_d = e_r*e_0
c_0 = 3e8

alpha = -w*np.sqrt(mu_0*e_d)*np.power(1+(sigma/(w*e_d))**2, 1/4)*np.sin(0.5*np.arctan(-sigma/(w*e_d)))
beta = w*np.sqrt(mu_0*e_d)*np.power(1+(sigma/(w*e_d))**2, 1/4)*np.cos(0.5*np.arctan(-sigma/(w*e_d)))
lbda = 2*np.pi/beta

k_c = beta - 1j*beta

z_10db = np.log(10)/alpha
print(z_10db)

f = np.logspace(-1, 20)
w = 2*np.pi*f
k_c = w*np.sqrt(mu_0*e_d)*np.power(1-1j*sigma/(w*e_d), 1/2)
delta = -1/np.imag(k_c)
plt.loglog(f, delta)
plt.show()
