import numpy as np
import scipy.optimize as sco
import matplotlib.pyplot as plt


Z_0 = 50
func = lambda R : 1-1/np.sqrt(2) - R**2/(R**2+R*Z_0+Z_0**2)
R_f = np.linspace(0, 300)
f = 1-1/np.sqrt(2) - R_f**2/(R_f**2+R_f*Z_0+Z_0**2)
plt.plot(R_f, f)
R = sco.fsolve(func, 50)
print(R/Z_0)
plt.show()
