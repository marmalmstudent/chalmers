import numpy as np
import matplotlib.pyplot as plt


load_data = np.transpose(np.loadtxt("data.txt", delimiter="\t"))
fig = plt.figure()
ax = plt.subplot(111)
ax.plot(load_data[0], load_data[1], label='$S_{11}$')
ax.plot(load_data[0], load_data[2], label='$S_{12}$')
ax.plot(load_data[0], load_data[3], label='$S_{21}$')
ax.plot(load_data[0], load_data[4], label='$S_{22}$')
ax.legend()
plt.show()
