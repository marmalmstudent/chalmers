import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
pi = np.pi
c_0 = 2.99792458e8
mu_0 = pi*4e-7
epsilon_0 = 1/(mu_0*np.power(c_0, 2))

kl = 160/180*pi
w = 2*pi*2e9
z_0 = 50
y_0 = 1/z_0
C = 0.5e-12
B = 1/(w*C)
e_r = 10  # assume
d = c_0/(w*np.sqrt(e_r))*kl
tl_mat1 = np.array([[np.cos(kl/2), 1j*z_0*np.sin(kl/2)],
                    [1j*y_0*np.sin(kl/2), np.cos(kl/2)]])
y_mat = np.array([[1, -1j*B],
                  [0, 1]])
tl_mat2 = np.array([[np.cos(kl/2), 1j*z_0*np.sin(kl/2)],
                    [1j*y_0*np.sin(kl/2), np.cos(kl/2)]])
abcd = np.dot(tl_mat1, np.dot(y_mat, tl_mat2))
print("\nABCD:\n", abcd)
print("Reciprocal: ", np.linalg.det(abcd) == 1)


def func_arr(kd):
    """
    func_arr(kd)

    Calculates (A+B)/2 (=cosh(gamma*d)) for a periodic structure
    of a transmission line, capacitor and a transmission line in series

    numpy.tensordot:
    Does dot product using tensors. This allows the matrix elements to be
    arrays for efficient computing. Two matrixes with shape (2,2,n) and (2,2,m)
    produces a matrixes with shape (2,m,2,n).
    numpy.diagonal:
    Since numpy.tensordot returns (2,m,2,n) matrix but we only want the
    parts where m=n, we take the diagonal with respect to axis 1 and 3. A
    matrix with shape (2,m,2,n), where m=n, produces a matrix with shape
    (2,2,n).
    numpy.transpose:
    Matrix axis 0 and 1 have switched place during during tensordot so we need
    to transpose it back to its original configuration. A matrix with shape
    (2,2,b) produces a matrix with shape (2,2,n)

    Parameter
    ---------
    kd : array_like
        This array is the electrical length of the line.

    Returns
    -------
    out : array
        From the resulting ABCD-matrix it returns (A+D)/2. This is equal to
        cosh(gamma*d) where gamma is the propagation constant.
    """
    tl_mat1 = np.array([[np.cos(kd/2), 1j*z_0*np.sin(kd/2)],
                        [1j*y_0*np.sin(kd/2), np.cos(kd/2)]])
    y_mat = np.array([[np.ones(np.size(kd)), -1j*np.sqrt(e_r)*d/(kd*c_0*C)],
                      [np.zeros(np.size(kd)), np.ones(np.size(kd))]])
    tl_mat2 = np.array([[np.cos(kd/2), 1j*z_0*np.sin(kd/2)],
                        [1j*y_0*np.sin(kd/2), np.cos(kd/2)]])
    tmp_mat = np.transpose(
        np.diagonal(
            np.tensordot(tl_mat1, y_mat, axes=(0, 1)),
            axis1=1, axis2=3),
        (1, 0, 2))
    tmp_mat2 = np.transpose(
        np.diagonal(
            np.tensordot(tmp_mat, tl_mat2, axes=(0, 1)),
            axis1=1, axis2=3),
        (1, 0, 2))
    out = np.real((tmp_mat2[0, 0]+tmp_mat2[1, 1])/2)
    return out


yl = np.linspace(-pi, pi, 101)
kd_1 = spo.fsolve(lambda kd: np.cos(yl) - func_arr(kd),
                  2.5*np.ones(np.size(yl)))
kd_2 = spo.fsolve(lambda kd: np.cos(yl) - func_arr(kd),
                  5*np.ones(np.size(yl)))
"""
kd_3 = spo.fsolve(lambda kd: np.cos(yl) - func_arr(kd),
                  (2*pi+1.5)*np.ones(np.size(yl)))
kd_4 = spo.fsolve(lambda kd: np.cos(yl) - func_arr(kd),
                  (4*pi-1)*np.ones(np.size(yl)))
"""
plt.plot(yl, kd_1, color="#000000")
plt.plot(yl, kd_2, color="#000000")
plt.plot([yl[0], yl[len(yl)-1]], [kd_1[0], kd_1[len(kd_1)-1]],
         linestyle="--", linewidth=1, color="#000000")
plt.plot([yl[0], yl[len(yl)-1]], [kd_2[0], kd_2[len(kd_2)-1]],
         linestyle="--", linewidth=1, color="#000000")
plt.plot([yl[0], yl[len(yl)-1]], [kd_1[np.argmin(kd_1)],
                                  kd_1[np.argmin(kd_1)]],
         linestyle="--", linewidth=1, color="#000000")
plt.plot([yl[0], yl[len(yl)-1]], [kd_2[np.argmax(kd_2)],
                                  kd_2[np.argmax(kd_2)]],
         linestyle="--", linewidth=1, color="#000000")
"""
plt.plot(yl, kd_3, color="#000000")
plt.plot(yl, kd_4, color="#000000")
"""
plt.xlabel("bl")
plt.ylabel("kl")
plt.axis([-pi, pi, 0, 7])
plt.show()
