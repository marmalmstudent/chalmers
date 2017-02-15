#! /usr/bin/python2
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sco
import numpy.matlib as npm
import os
pi = np.pi


def plotOctave(X, Y, Z):
    np.savetxt('X.txt', X)
    np.savetxt('Y.txt', Y)
    np.savetxt('Z.txt', Z)
    with open("OctavePlotRun.m", "w") as text_file:
        print("[hFig, hSurf] = OctavePlot(load('X.txt'), \"X\",",
              "load('Y.txt'), \"Y\", load('Z.txt'), \"Z\", 1);\n",
              "while(ishandle(hFig)) # sleep until figure is closed\n",
              "sleep(0.5);\n",
              "end\n",
              "exit(0) # close octave and free up some memory\n",
              file=text_file)
        # print("Purchase Amount: {}".format(TotalAmount), file=text_file)
    os.system('octave --no-window-system --silent --persist OctavePlotRun.m &')


W = 5e-3
# f = 30e9
# f = 44.4e9
f = 50e9
w = 2*np.pi*f
mu_0 = 4*np.pi*1e-7
epsilon_0 = 8.85e-12
epsilon_r = 3
func = lambda K_d : np.sqrt(
    w**2*mu_0*epsilon_0*(epsilon_r-1)-K_d)*np.tan(K_d*W/2) - K_d
fig = plt.figure()
margin = 1e1
bg_clr = '#ffffff'
line_clr = '#000000'
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel('')
ax.set_ylabel('')
ax.grid(b=True, color=line_clr)
ax.set_axis_bgcolor(bg_clr)
K_d_arr = np.linspace(0, pi/W-margin, 501)
np.savetxt('K_d_arr1.txt', K_d_arr)
np.savetxt('func1.txt', func(K_d_arr))
ax.plot(K_d_arr, func(K_d_arr), color=line_clr)
K_d_arr = np.linspace(pi/W+margin, 3*pi/W-margin, 1001)
np.savetxt('K_d_arr2.txt', K_d_arr)
np.savetxt('func2.txt', func(K_d_arr))
ax.plot(K_d_arr, func(K_d_arr), color=line_clr)
ax.axis([0, 3*pi/W-margin, func(-pi/W+margin), func(3*pi/W-margin)])

K_d_initial_guess = 1700
K_d = sco.fsolve(func, K_d_initial_guess)
np.savetxt('K_d.txt', K_d)
np.savetxt('func.txt', func(K_d))
ax.plot(K_d, func(K_d), linestyle='', linewidth=5,
        marker='D', markersize=5, color='#ff0000')
K_a = np.sqrt(w**2*mu_0*epsilon_0*(epsilon_r-1)-K_d)
gamma_d_sqr = K_d**2 - w**2*mu_0*epsilon_0*epsilon_r
if (gamma_d_sqr >= 0):
    gamma_d = np.sqrt(gamma_d_sqr)
else:
    gamma_d = 1j*np.sqrt(-gamma_d_sqr)
gamma_0_sqr = K_d**2 - w**2*mu_0*epsilon_0
if (gamma_0_sqr >= 0):
    gamma_0 = np.sqrt(gamma_0_sqr)
else:
    gamma_0 = 1j*np.sqrt(-gamma_0_sqr)
f_c = K_d/(2*np.pi*np.sqrt(mu_0*epsilon_0*epsilon_r))
print(K_d, K_a, f_c)
# print(K_d, K_a, gamma_d, gamma_0, f_c)
E_0y = 1e9

z = np.linspace(0, 5e-3, 101)
x_01 = np.linspace(-3*W/2, -3*W/2 + W*99/100, 100)
x_d = np.linspace(-W/2, W/2, 101)
x_02 = np.linspace(W/2+W/100, W/2 + W, 100)
x = np.array([])
x = np.append(x, x_01)
x = np.append(x, x_d)
x = np.append(x, x_02)
Z = npm.repmat(z, len(x), 1)
X = npm.repmat(x, len(z), 1)

H_z01 = -1j*K_d**2*E_0y/(
    w*mu_0*K_d)*np.sin(K_d*W/2)*np.exp(-K_a*(np.abs(x_01)-W/2))
H_zd = -1j*K_d**2*E_0y/(w*mu_0*K_d)*np.sin(K_d*x_d)
H_z02 = -1j*K_d**2*E_0y/(
    w*mu_0*K_d)*np.sin(K_d*W/2)*np.exp(-K_a*(np.abs(x_02)-W/2))
H_z = np.array([])
H_z = np.append(H_z, H_z01)
H_z = np.append(H_z, H_zd)
H_z = np.append(H_z, H_z02)
H_Z = np.abs(np.real(H_z*np.exp(-gamma_d*Z).T))

H_x01 = -1j*np.sign(x_01)*gamma_d*E_0y/(
    w*mu_0)*np.cos(K_d*W/2)*np.exp(-K_a*(np.abs(x_01)-W/2))
H_xd = -1j*gamma_d*E_0y/(w*mu_0)*np.cos(K_d*x_d)
H_x02 = -1j*np.sign(x_02)*gamma_d*E_0y/(
    w*mu_0)*np.cos(K_d*W/2)*np.exp(-K_a*(np.abs(x_02)-W/2))
H_x = np.array([])
H_x = np.append(H_x, H_x01)
H_x = np.append(H_x, H_xd)
H_x = np.append(H_x, H_x02)
H_X = np.abs(np.real(H_x*np.exp(-gamma_d*Z).T))

E_y01 = np.sign(
    x_01)*E_0y*np.cos(K_d*W/2)*np.exp(-K_a*(np.abs(x_01)-W/2))
E_yd = E_0y*np.cos(K_d*x_d)
E_y02 = np.sign(
    x_01)*E_0y*np.cos(K_d*W/2)*np.exp(-K_a*(np.abs(x_02)-W/2))
E_y = np.array([])
E_y = np.append(E_y, E_y01)
E_y = np.append(E_y, E_yd)
E_y = np.append(E_y, E_y02)
E_Y = np.abs(np.real(E_y*np.exp(-gamma_d*Z).T))

H = np.abs(np.real(H_x*np.exp(-gamma_d*Z).T)-np.real(H_z*np.exp(-gamma_d*Z).T))
H = H*(H < 1e5) + 1e5*(H >= 1e5)
print(K_d, gamma_d)

"""
fig1 = plt.figure()
ax1 = fig1.gca(projection='3d')
ax1.plot_surface(X, Z.T, H,
                 rstride=1, cstride=1, cmap=mpl.cm.jet,
                 linewidth=0, antialiased=False)
ax1.set_ylabel('z')
ax1.set_xlabel('x')
ax1.set_zlabel('abs(real(H_x)-real(H_z))')
# ax1.view_init(elev=30, azim=-50)
ax1.view_init(elev=30, azim=66)
ax1.dist = 8
"""
# print(np.size(X_01, axis=0), np.size(Z_01, axis=0), np.size(H_z01*np.exp(-gamma_d*Z_01).T, axis=0))
# print(np.size(X_01, axis=1), np.size(Z_01, axis=1), np.size(H_z01*np.exp(-gamma_d*Z_01).T, axis=1))
"""
ax1 = fig1.add_subplot(2, 2, 1)
ax1.plot(x, np.abs(np.real(H_z)), color='#ff0000')

ax2 = fig1.add_subplot(2, 2, 2)
ax2.plot(x, np.abs(np.real(H_x)), color='#00ff00')

ax3 = fig1.add_subplot(2, 2, 3)
ax3.plot(x, np.abs(np.real(E_y)), color='#0000ff')
"""

plotOctave(X, Z.T, E_Y)
#plt.show()
