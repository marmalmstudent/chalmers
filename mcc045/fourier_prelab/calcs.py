import numpy as np
import numpy.fft as npf
import scipy.integrate as spi
import scipy.optimize as spo
import matplotlib.pyplot as plt
"""
z = np.linspace(0, 6.05)
lbda = 632.8e-9
w_z = 25/16*1e-3
print(np.sqrt(w_z**4/4-(z*lbda/np.pi)**2), w_z**2/2)
w_0 = np.sqrt(w_z**2/2+np.sqrt(w_z**4/4-(z*lbda/np.pi)**2))
plt.plot(z, w_0)
# print(w_0)
# plt.axis([0, 1, 0, w_z])
plt.show()
"""

print(spo.fsolve(lambda x: 1/(400e-3-x)+1/(400e-3-x)-1/200e-3, 399e-3))


"""
def integrand(a):
    return spi.quad(lambda x: np.exp(-x**2), -a, a)


#x=np.linspace(-5,5,201)
#y=np.exp(-x**2))
bounds = spo.fsolve(lambda a: integrand(a)-np.exp(-2), 1)
arr = np.linspace(-bounds, bounds, 1001)
func = np.exp(-arr**2)
freqs = npf.fft(func)
plt.plot(np.abs(freqs))
plt.show()
"""
"""
plt.plot(x,y)
plt.show()
"""
