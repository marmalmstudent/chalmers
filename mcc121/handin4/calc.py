import numpy as np


S = np.array(
    [[0.2*np.exp(1j*np.pi/2), 0.5*np.exp(-1j*np.pi/4), 0.5*np.exp(1j*np.pi/4), 0],
     [0.5*np.exp(-1j*np.pi/4), 0, 0, 0.5*np.exp(1j*np.pi/4)],
     [0.5*np.exp(-1j*np.pi/4), 0, 0, 0.5*np.exp(-1j*np.pi/4)],
     [0, 0.5*np.exp(-1j*np.pi/4), 0.5*np.exp(-1j*np.pi/4), 0]])
result1 = np.dot(S, np.conj(S.T))
print(result1)
