import numpy as np
import matplotlib.pyplot as plt

mat = np.loadtxt('DOE.txt', comments='#',
                 delimiter=None, converters=None, skiprows=5,
                 usecols=None, unpack=False, ndmin=0)
print(np.size(mat))
fig = plt.figure()
ax = fig.add_subplot(111)
ax.matshow(mat, cmap=plt.cm.gray)
plt.show()
