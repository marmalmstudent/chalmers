import numpy as np

G = np.power(10, np.array([-0.5, 10, -4, -0.5, 20])/10)
F = np.power(10, np.array([0.5, 2, 5, 0.5, 3])/10)
G = np.array([-0.5, 10, -4, -0.5, 20])
F = np.array([0.5, 2, 5, 0.5, 3])
p = 0
print(np.power(10, (0 - G[3] - G[1] - G[2] - G[0] - G[4] - G[0])/10))
"""
F_sys_rx = F[0] + (F[1]-1)/(G[0]) + (F[2]-1)/(G[0]*G[1]) + (F[3]-1)/(G[0]*G[1]*G[2]) + (F[4]-1)/(G[0]*G[1]*G[2]*G[3])

F_sys_tx = F[3] + (F[1]-1)/(G[3]) + (F[2]-1)/(G[3]*G[1]) + (F[0]-1)/(G[3]*G[1]*G[2]) + (F[4]-1)/(G[3]*G[1]*G[2]*G[0]) + (F[0]-1)/(G[3]*G[1]*G[2]*G[0]*G[4])
F_sys_rx_db=10*np.log10(F_sys_rx)
F_sys_tx_db=10*np.log10(F_sys_tx)
print("F_rx: ", F_sys_rx, "F_rx[dB]: ", F_sys_rx_db)
print("F_tx", F_sys_tx, "F_tx[dB]: ", F_sys_tx_db)
"""
