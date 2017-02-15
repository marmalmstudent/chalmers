import numpy as np
import matplotlib.pyplot as plt

Tb2 = 3.9
Tb1 = 17.7
Tbg = 2.7
Tex = np.linspace(25, 25.5, 2001)
f1 = ((Tb2-Tbg)/(Tbg-Tex)+1)**20
f2 = ((Tb1-Tbg)/(Tbg-Tex)+1)
print((Tb1-Tbg)/(Tbg-25.26))
print((Tb2-Tbg)/(Tbg-25.26))
plt.plot(Tex, f1-f2)
plt.plot(Tex, np.zeros(np.size(Tex, axis=0)))
# plt.yscale('log')
plt.show()
