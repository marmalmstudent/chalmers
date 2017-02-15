import numpy as np
import scipy.optimize as spo

def func(x):
    x1, x2, lbda1, lbda2 = x
    f1 = 1 + 2*lbda1*(x1-1) - 2*lbda2*x1
    f2 = 2*lbda1*(x2+2)-2*lbda2*x2
    f3 = lbda1*((x1-1)**2 + (x2+2)**2 - 16)
    f4 = lbda2*(13-x1**2-x2**2)
    return f1, f2, f3, f4


x_1 = np.linspace(-5, 6, 10)
x_2 = np.linspace(-7, 4, 10)
for i in range(0, 10):
    for j in range(0, 10):
        print(spo.fsolve(lambda x: func(x), [i, j, 1, 1]))
