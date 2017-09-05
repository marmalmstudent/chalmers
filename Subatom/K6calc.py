import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ~~~~~angle
d1 = 2.5  # cm
d2 = 5.5  # cm
d3 = 8.5  # cm
theta1 = 1
theta2 = (d1/d2)**2
theta3 = (d1/d3)**2
# ~~~~~angle
# ~~~~~
x = np.array([35, 55, 75, 1*60+35, 1*60+55, 2*60+15, 2*60+35, 2*60+55, 3*60+15, 3*60+35, 3*60+55, 4*60+20, 4*60+50, 5*60+20, 5*60+50, 6*60+20, 6*60+50, 7*60+20, 7*60+50, 8*60+20, 8*60+50, 9*60+20, 9*60+50, 10*60+35, 11*60+35, 12*60+35])  # time
y = np.array([1107, 895, 707, 514, 397, 328, 270, 251, 210, 195, 163, 307/2.0, 207/2.0, 245/2.0, 196/2.0, 175/2.0, 133/2.0, 146/2.0, 119/2.0, 102/2.0, 103/2.0, 73/2.0, 62/2.0, 132/5.0, 90/5.0, 74/5.0])*6-14  # counts/min minus background
eq1 = np.polyfit(x[11:len(x)], np.log(y[11:len(x)]), 1)  # long-lived isotopes qeuation (lambda*x1+x0)
lbda1 = eq1[0]  # lambda
A1 = np.exp(eq1[1])
y_ll = A1*np.exp(lbda1*x)
y_sl_temp = (y-y_ll)[0:6]
print np.log(2.0)/lbda1, A1
eq2 = np.polyfit(x[0:6], np.log(y_sl_temp[0:6]), 1)  # long-lived isotopes qeuation (lambda*x1+x0)
lbda2 = eq2[0]
A2 = np.exp(eq2[1])
y_sl = A2*np.exp(lbda2*x)
y_tot = y_ll+y_sl
print np.log(2.0)/lbda2, A2
eq3 = np.polyfit(x[0:6], np.log(y[0:6]), 1)  # long-lived isotopes qeuation (lambda*x1+x0)
lbda3 = eq3[0]
A3 = np.exp(eq3[1])
y_llsl = A3*np.exp(lbda3*x)
print lbda3, A3

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
line = ax.plot(x, y, label='Measured total activity', linestyle='', marker='o', color='k')  # measured data
line1 = ax.plot(x, y_ll, label='Long-lived isotope acitivity', linestyle='-', marker=None, color='r')
line2 = ax.plot(x[0:8], y_llsl[0:8], label='Long-lived and short-lived activity', linestyle='-', marker=None, color='g')
line21 = ax.plot(x[0:12], y_sl[0:12], label='Short-lived isotope acitivity', linestyle='--', marker=None, color='g')
line3 = ax.plot(x, y_tot, label='Total activity', linestyle='-', marker=None, color='b')
plt.legend(numpoints=1, loc=1)  #  handles=[line, line1, line2, line21, line3]
ax.set_yscale('log')
plt.axis([0, 800, 50, 1e4])
plt.ylabel('Number of active Ag')
plt.xlabel('Time, [s]')
plt.show()
