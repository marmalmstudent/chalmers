import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sco


def smith_chart_circles_real():
    """
    r = np.array([[0.1], [0.2], [0.3], [0.5], [0.5], [0.6], [0.7], [0.8],
                  [0.9], [1.0], [1.2], [1.4], [1.6], [1.8], [2.0], [3.0],
                  [4.0], [5.0], [10.0], [20.0], [50.0]])
    """
    r = np.array([[0.2], [0.5], [1], [2], [5]])
    x = []
    y = []
    for i in range(0, np.size(r, axis=0)):
        theta = 2*np.pi - np.linspace(0, 2*np.pi, int(1000/(1+r[i])))
        x_temp = r[i]/(1+r[i]) + 1/(1+r[i])*np.cos(theta)
        y_temp = 1/(1+r[i])*np.sin(theta)
        x.append(x_temp)
        y.append(y_temp)
    return x, y


def smith_chart_circles_imag():
    """
    X = np.array([[-50.0], [-20.0], [-10.0], [-5.0], [-4.0], [-3.0], [-2.0],
                  [-1.8], [-1.6], [-1.4], [-1.2], [-1.0], [-0.9], [-0.8],
                  [-0.7], [-0.6], [-0.5], [-0.4], [-0.3], [-0.2], [-0.1],
                  [0.1], [0.2], [0.3], [0.4], [0.5], [0.6], [0.7], [0.8],
                  [0.9], [1.0], [1.2], [1.4], [1.6], [1.8], [2.0], [3.0],
                  [4.0], [5.0], [10.0], [20.0], [50.0]])
    """
    X_1 = np.array([[0.2], [0.5], [1], [2], [5]])
    X_2 = np.array([[-5], [-2], [-1], [-0.5], [-0.2]])
    x = []
    y = []
    for i in range(0, np.size(X_1, axis=0)):
        theta_1 = np.linspace(np.pi/2, 3*np.pi/2, int(1000/np.abs(X_1[i])))
        x_temp_1 = 1 + np.cos(theta_1)/X_1[i]
        y_temp_1 = 1/X_1[i] + np.sin(theta_1)/X_1[i]
        x_1 = np.array([])
        y_1 = np.array([])
        for j in range(0, np.size(x_temp_1)-1):
            if (x_temp_1[j]**2 + y_temp_1[j]**2 <= 1):
                x_1 = np.append(x_1, x_temp_1[j])
                y_1 = np.append(y_1, y_temp_1[j])
        x.append(x_1)
        y.append(y_1)
    for i in range(0, np.size(X_2, axis=0)):
        theta_2 = np.linspace(3*np.pi/2, 5*np.pi/2, int(1000/np.abs(X_2[i])))
        x_temp_2 = 1 + np.cos(theta_2)/X_2[i]
        y_temp_2 = 1/X_2[i] + np.sin(theta_2)/X_2[i]
        x_2 = np.array([])
        y_2 = np.array([])
        for j in range(0, np.size(x_temp_2)-1):
            if (x_temp_2[j]**2 + y_temp_2[j]**2 <= 1):
                x_2 = np.append(x_2, x_temp_2[j])
                y_2 = np.append(y_2, y_temp_2[j])
        x.append(x_2)
        y.append(y_2)
    return x, y


imp_clr = '#ff0000'
adm_clr = '#0000ff'
normal_clr = '#000000'
intermediate_imp_style = '--'
imp_change_style = ':'
final_imp_style = '-'


fig1 = plt.figure()
ax1 = fig1.add_subplot(111, title='Inerted-F antenna matching')
smith_chart_linewidth = 0.3
a, b = smith_chart_circles_real()
for i in range(0, np.size(a, axis=0)):
    ax1.plot(np.array(a[i]), np.array(b[i]), linewidth=smith_chart_linewidth,
             linestyle='-', marker='', markersize=1, color=imp_clr)
    ax1.plot(-np.array(a[i]), -np.array(b[i]), linewidth=smith_chart_linewidth,
             linestyle='-', marker='', markersize=1, color=adm_clr)

c, d = smith_chart_circles_imag()
for j in range(0, np.size(c, axis=0)):
    ax1.plot(np.array(c[j]), np.array(d[j]), linewidth=smith_chart_linewidth,
             linestyle='-', marker='', markersize=1, color=imp_clr)
    ax1.plot(-np.array(c[j]), np.array(d[j]), linewidth=smith_chart_linewidth,
             linestyle='-', marker='', markersize=1, color=adm_clr)

ax1.plot(np.linspace(-1, 1, 2), np.zeros(2), linewidth=smith_chart_linewidth,
         linestyle='-', marker='', markersize=1, color=imp_clr)
ax1.plot(np.cos(np.linspace(0, 2*np.pi, 1000)),
         np.sin(np.linspace(0, 2*np.pi, 1000)),
         linewidth=smith_chart_linewidth, linestyle='-', marker='',
         markersize=1, color=imp_clr)

Z_0 = 50  # system impedance

# Load impedance
Z_L = 60 - 100j
F_0 = (Z_L-Z_0)/(Z_L+Z_0)
ax1.plot(np.real(F_0), np.imag(F_0), linewidth=2,
         linestyle=final_imp_style, marker='o',
         markersize='2', color=normal_clr)

# Transmission line
func = lambda v: np.real(1/(Z_0*(1+np.abs(F_0)*np.exp(1j*v))/(1-np.abs(F_0)*np.exp(1j*v)))) - 1/Z_0
v = sco.fsolve(func, -np.imag(F_0))
Z_1 = Z_0*(1+np.abs(F_0)*np.exp(1j*v))/(1-np.abs(F_0)*np.exp(1j*v))

F_1 = (Z_1-Z_0)/(Z_1+Z_0)
ax1.plot(np.real(F_1), np.imag(F_1), linewidth=2,
         linestyle=final_imp_style, marker='o',
         markersize='2', color=normal_clr)
arg_F_0 = (2*np.pi + np.angle(F_0)) % (2*np.pi)
arg_F_1 = (2*np.pi + np.angle(F_1)) % (2*np.pi)
v_tmp = np.linspace(np.min((arg_F_0, arg_F_1)),
                    np.max((arg_F_0, arg_F_1)))
F_1_1 = np.abs(F_0)*np.exp(1j*v_tmp)
ax1.plot(np.real(F_1_1), np.imag(F_1_1), linewidth=1,
         linestyle=intermediate_imp_style, marker='',
         markersize='2', color=normal_clr)
TL_bl = np.max((arg_F_0, arg_F_1)) -\
     np.min((arg_F_0, arg_F_1))
print('TL length: ', TL_bl/(2*np.pi)*(np.pi)/(2*np.pi), ' lambda')

# Opend stub
Z_open = sys.float_info.max / 1e10
func1 = lambda v: np.imag(1/(Z_0*(Z_open+1j*Z_0*np.tan(v))/(Z_0+1j*Z_open*np.tan(v)))) + np.imag(1/Z_1)
stub_bl = sco.fsolve(func1, np.pi/3)
Z_stub = Z_0*(Z_open+1j*Z_0*np.tan(stub_bl))/(Z_0+1j*Z_open*np.tan(stub_bl))
F_stub = (Z_stub - Z_0)/(Z_stub + Z_0)
ax1.plot(np.real(F_stub), np.imag(F_stub), linewidth=1,
         linestyle=intermediate_imp_style, marker='o',
         markersize='2', color=normal_clr)
bl = np.linspace(0, stub_bl)
Z_stub_1 = Z_0*(Z_open+1j*Z_0*np.tan(bl))/(Z_0+1j*Z_open*np.tan(bl))
F_stub_1 = (Z_stub_1 - Z_0)/(Z_stub_1 + Z_0)
ax1.plot(np.real(F_stub_1), np.imag(F_stub_1), linewidth=1,
         linestyle=intermediate_imp_style, marker='',
         markersize='2', color=normal_clr)

Z_2 = 1/(1/Z_1 + 1/(Z_0*(Z_open+1j*Z_0*np.tan(stub_bl))/(Z_0+1j*Z_open*np.tan(stub_bl))))
F_2 = (Z_2 - Z_0)/(Z_2 + Z_0)
ax1.plot(np.real(F_2), np.imag(F_2), linewidth=1,
         linestyle=intermediate_imp_style, marker='o',
         markersize='5', color=normal_clr)
Z_2_1 = 1/(1/Z_1 + 1/(Z_0*(Z_open+1j*Z_0*np.tan(bl))/(Z_0+1j*Z_open*np.tan(bl))))
F_2_1 = (Z_2_1 - Z_0)/(Z_2_1 + Z_0)
ax1.plot(np.real(F_2_1), np.imag(F_2_1), linewidth=1,
         linestyle=intermediate_imp_style, marker='',
         markersize='2', color=normal_clr)
arg_F_stub = (2*np.pi + np.angle(F_stub)) % (2*np.pi)
S_bl = 2*np.pi - arg_F_stub
print('TL length: ', S_bl/(2*np.pi)*(np.pi)/(2*np.pi), ' lambda')


"""
# transmission line part
Z_1 = Z_0 + R  # at lambda/4 end looking at R
F_1 = (Z_1-Z_0)/(Z_1+Z_0)

ax1.plot((np.real(F_0), np.real(F_1)), (np.imag(F_0), np.imag(F_1)),
         linewidth=1, linestyle=imp_change_style, color=normal_clr)
ax1.plot(np.real(F_1), np.imag(F_1), linewidth=2,
         linestyle=final_imp_style, marker='o',
         markersize='2', color=normal_clr)

bl = np.linspace(0, np.pi/2)
# at R looking at lambda/4 end
Z_2_0 = Z_1
F_2_0 = (Z_2_0 - Z_0)/(Z_2_0 + Z_0)
ax1.plot(np.real(F_2_0), np.imag(F_2_0), linewidth=1,
         linestyle=intermediate_imp_style, color=normal_clr)
Z_2_1 = Z_0*(Z_1 +
             1j*Z_0*np.tan(bl))/(Z_0 + 1j*Z_1*np.tan(bl))
F_2_1 = (Z_2_1-Z_0)/(Z_2_1+Z_0)
ax1.plot(np.real(F_2_1), np.imag(F_2_1), linewidth=1,
         linestyle=imp_change_style, color=normal_clr)

bl = np.pi/2
Z_2 = Z_0*(Z_1 + 1j*Z_0*np.tan(bl))/(Z_0 + 1j*Z_1*np.tan(bl))
F_2_3 = (Z_2 - Z_0)/(Z_2 + Z_0)
F_2 = (Z_2-R)/(Z_2+R)
ax1.plot(np.real(F_2_3), np.imag(F_2_3), linewidth=1,
         linestyle=intermediate_imp_style, marker='o',
         markersize='2', color=normal_clr)
ax1.plot((np.real(F_2_3), np.real(F_2)), (np.imag(F_2_3), np.imag(F_2)),
         linewidth=1, linestyle=imp_change_style, color=normal_clr)
ax1.plot(np.real(F_2), np.imag(F_2), linewidth=2,
         linestyle=final_imp_style, marker='o',
         markersize='2', color=normal_clr)
Z_3 = R + Z_2  # at Z_c lookint at R
F_3 = (Z_0-Z_3)/(Z_0+Z_3)
ax1.plot((np.real(F_2), np.real(F_3)), (np.imag(F_2), np.imag(F_3)),
         linewidth=1, linestyle=imp_change_style, color=normal_clr)
ax1.plot(np.real(F_3), np.imag(F_3), linewidth=2,
         linestyle=final_imp_style, marker='o',
         markersize='5', color=normal_clr)


Gamma = F_3
RL = -20*np.log10(np.abs(Gamma))
VSWR = (1+np.abs(Gamma))/(1-np.abs(Gamma))
print(Gamma, RL, VSWR)
"""
ax1.set_aspect('equal')
ax1.axis([-1.1, 1.1, -1.1, 1.1])
plt.show()
